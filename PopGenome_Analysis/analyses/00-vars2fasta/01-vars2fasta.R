#----------------------------------------------------------------------------------------------------
# Purpose of this script is to filter samples with too many heterozygous calls
# and to force het calls to major allele to convert to fasta
#----------------------------------------------------------------------------------------------------
library(seqinr)
library(vcfR)
library(vcfRmanip)
library(tidyverse)
source("R/vcf_utils.R")
set.seed(48)
#........................................................
# Read in metadata
#........................................................
smpls <- readxl::read_excel("~/Documents/GitHub/VivID_Seq/scrape_pubseqs/vivid_seq_public_NGS.xlsx") %>%
  dplyr::select(c("iid", "acc", "country", "vividregion"))

drc <- data.frame(iid = c("D9U3K", "O6Y4K", "Q8J6O"),
                  acc = c("D9U3K", "O6Y4K", "Q8J6O"),
                  country = "CD",
                  vividregion = "AF")

smpls <- rbind.data.frame(smpls, drc) %>%
  dplyr::mutate(country = ifelse(vividregion == "NHA", "NHA", country)) %>%
  dplyr::filter(vividregion != "Lab")

#........................................................
# Read in fastas
#........................................................
pvp01 <- seqinr::read.fasta("~/Documents/MountPoints/mountIDEEL/resources/genomes/Pvivax/genomes/PvP01.fasta",
                            forceDNAtolower = F)
pvp01 <- pvp01$PvP01_MIT_v1

#........................................................
# Read in vcf
#........................................................
vcf <- vcfR::read.vcfR(file = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/vcfs_hardfilt_variants/passed.joint.vcf.gz")
vcf.snp <- vcfR::extract.indels(vcf, return.indels = F)
vcf.snp

#........................................................
# Find samples with a high level of heterozygosity
# which mimicks COI > 1 but is really heteroplasmy
# or sequencing error
#........................................................
smplhet <- apply(vcfR::is.het(vcfR::extract.gt(vcf.snp)), 2, mean)
smplhet <- tibble::tibble(smpls = names(smplhet),
                          hetprop = smplhet)
smplhetpass <- smplhet$smpls[smplhet$hetprop < 0.2] # less than 20%

vcf.snp <- vcfRmanip::select_samples(vcfRobject = vcf.snp,
                                     smplvctr = smplhetpass)

#........................................................
# Remove loci that have deleted allele alternative
# which is encoded in vcf as a *
#........................................................
alts <- vcfR::getALT(vcf.snp)
asterickfilt <- grepl(pattern = "\\*", x = alts)
asterickfilt <- data.frame(seqname = vcfR::getCHROM(vcf.snp),
                           start = vcfR::getPOS(vcf.snp),
                           end = vcfR::getPOS(vcf.snp),
                           geneid = 1:length(vcfR::getPOS(vcf.snp))) %>%
  dplyr::filter(asterickfilt)

vcf.snp <- vcfRmanip::vcffilter_ChromPos(vcfRobject = vcf.snp,
                                         chromposbed = asterickfilt)


#........................................................
# Impute Missing Calls
#........................................................
# split by country for imputation
psdsmpls <- tibble::tibble(acc = smplhetpass)
psdsmpls <- dplyr::inner_join(psdsmpls, smpls) # take into account samples we just dropped

# NOTE, for imputation purpose LA (Laos is going to be combined into Thailand)
# NOTE, for imputation purpose LK (Sri Lanka is going to be combined into India)
psdsmpls$country[psdsmpls$acc == "ERS347479"] <- "TH" # instance of LA
psdsmpls$country[psdsmpls$acc == "ERS040109"] <- "IN" # instance of LK

smpls.split <- split(psdsmpls$acc, psdsmpls$country)
mtdna.list <- lapply(smpls.split, function(x){
  ret <- vcfRmanip::select_samples(vcfRobject = vcf.snp,
                                   smplvctr = x)
  return(ret) })


mtdna.list.imp <- lapply(mtdna.list, vcfR_impute_gt)

#...........................
# Clean up imputation
#...........................
vcf.snp.imp <- cbindvcflist(mtdna.list.imp)

#........................................................
# Clean up het calls
#........................................................
vcf.snp.hom <- vcfR_force_hets_to_homo(vcf.snp.imp)
vcf.snp.hom <- vcfR_filter_unique_sites(vcf.snp.hom)

#........................................................
# Remove Private Alleles within a country
#........................................................
vcf.snp.hom.long <- vcf.snp.hom %>%
  vcfR::extract_gt_tidy(., format_fields = "GT") %>%
  dplyr::group_by(Indiv) %>%
  dplyr::mutate(POS = vcfR::getPOS(vcf.snp.hom)) %>%
  dplyr::rename(acc = Indiv) %>%
  dplyr::left_join(., smpls, by = "acc") %>%
  dplyr::ungroup()

countryn <- vcf.snp.hom.long %>%
  dplyr::select(c("acc", "country")) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(countryn = n())

privateloci.withincntry <- vcf.snp.hom.long %>%
  dplyr::left_join(., countryn, by = "country") %>%
  dplyr::mutate(I = 1,
                limit = 0.1 * (countryn)/(countryn - 1), # note correction for small sample size
                limit2 = 1/countryn
  ) %>%
  dplyr::group_by(country, POS, gt_GT) %>%
  dplyr::summarise(
    locicount = sum(I)/mean(countryn), # countryn all same
    limit = mean(limit), # limit all same
    limit2 = mean(limit2) # limit all same
  ) %>%
  dplyr::mutate(rm = ifelse(locicount <= limit & !is.na(gt_GT) & limit != Inf, T, F), # note catch for single sample
                rm = ifelse(locicount == limit2 & !is.na(gt_GT) & limit2 != 1, T, rm) # note catch for single sample
  ) %>%
  dplyr::select(c("country", "POS", "gt_GT", "rm"))


vcf.snp.hom.long.nopriv <- vcf.snp.hom.long %>%
  dplyr::left_join(., y = privateloci.withincntry, by = c("country", "POS", "gt_GT")) %>%
  dplyr::filter(rm != T) %>%
  dplyr::select(-c("rm"))


#........................................................
# Final Object
#........................................................
saveRDS(object = vcf.snp.hom.long.nopriv,
        file = "data/derived_data/final_vcf_snps.RDS")

#........................................................
# Liftover to sequence call
#........................................................
mksngbp <- function(x){
  x <- stringr::str_split_fixed(x, "/", n=2)[,1]
  return(x)
}

vcf.snp.hom.long.nopriv$gt_GT_alleles_sbp <- sapply(vcf.snp.hom.long.nopriv$gt_GT_alleles,
                                                    mksngbp)


#........................................................
# Make new fastas
#........................................................
seqlist <- vcf.snp.hom.long.nopriv %>%
  dplyr::select(c("acc", "POS", "gt_GT_alleles_sbp"))
seqlist <- split(seqlist, factor(seqlist$acc))


mkmut <- function(x){
  mut <- pvp01
  mut[x$POS] <- x$gt_GT_alleles_sbp
  attr(mut, "name") <- unique(x$acc)
  return(mut)

}

seqlist <- lapply(seqlist, mkmut)



#........................................................
# Read in and append Pcynomolgi
#........................................................
Pcynomolgi <- seqinr::read.fasta("PATH",
                                 forceDNAtolower = F)
finalfastas <- seqlist
finalfastas[[length(finalfastas) + 1 ]] <- Pcynomolgi[[1]]



#........................................................
# From the Gilabert Manuscript, the Pv-like Clade 2
# samples ERS333076 and ERS352725 have an extra order
# of magnitude more diversity than the other global isolates
# going to go ahead and exclude them
#........................................................
finalfastas <- finalfastas[ ! names(finalfastas) %in% c("ERS333076", "ERS352725") ]

# now have to write out final sample list
outsmpls <- smplhetpass[! smplhetpass %in% c("ERS333076", "ERS352725") ]
saveRDS(outsmpls, file = "data/derived_data/final_smpl_list.RDS")



#........................................................
# Write Out Final Fasta
#........................................................
seqnames <- c(names(finalfastas))
seqnames[length(seqnames)] <- "Pcynomolgi"
dir.create(path = "data/fasta/", recursive = T)
seqinr::write.fasta(sequences = finalfastas,
                    names = seqnames, file.out = "data/fasta/mtdna_anc.fa")



