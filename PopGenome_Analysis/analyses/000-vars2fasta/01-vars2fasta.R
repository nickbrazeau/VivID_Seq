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
smpls <- readxl::read_excel("../scrape_pubseqs/vivid_seq_public_NGS.xlsx") %>%
  dplyr::select(c("iid", "acc", "country", "vividregion"))

drc <- data.frame(iid = c("D9U3K", "O6Y4K", "Q8J6O"),
                  acc = c("D9U3K", "O6Y4K", "Q8J6O"),
                  country = "CD",
                  vividregion = "AF")

smpls <- rbind.data.frame(smpls, drc) %>%
  dplyr::mutate(country = ifelse(vividregion == "NHA", "NHA", country))

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
vcf.snp <- vcf.snp[vcfR::is.biallelic(vcf.snp), ]

#........................................................
# Find samples that could not be resolved by joint
# genotyping and had more than 5% missing
#........................................................
smpl_missinglocivcf <- vcfRmanip::calc_loci_missingness_by_smpl(vcf.snp)
smpl_missingpass <- smpl_missinglocivcf$sample[smpl_missinglocivcf$missprop < 0.05]
vcf.snp <- vcfRmanip::select_samples(vcfRobject = vcf.snp,
                                     smplvctr = smpl_missingpass)

#........................................................
# Find samples with a high level of heterozygosity
# which mimicks COI > 1 but is really heteroplasmy
# or sequencing error
#........................................................
smplhet <- apply(vcfR::is.het(vcfR::extract.gt(vcf.snp)), 2, mean)
smplhet <- tibble::tibble(smpls = names(smplhet),
                          hetprop = smplhet)
smplhetpass <- smplhet$smpls[smplhet$hetprop < 0.15] # less than 15%

vcf.snp <- vcfRmanip::select_samples(vcfRobject = vcf.snp,
                                     smplvctr = smplhetpass)


#........................................................
# Final Objects
#........................................................
dir.create("data/derived_data/")
saveRDS(object = vcf.snp,
        file = "data/derived_data/final_vcf_snps.RDS")

#............................................................
# coerce to single call (i.e. monoclonal)
#...........................................................
vcf.snp.mncnl <- vcfR_force_hets_to_homo(vcfRobj = vcf.snp)
# ~0 missing is correct based on extract gt funcion
sum(is.na(extract.gt(vcf.snp)))/length(extract.gt(vcf.snp))

# now make long
vcf.snp.mncnl.long <- vcf.snp.mncnl %>%
  vcfR::extract_gt_tidy(., format_fields = "GT") %>%
  dplyr::group_by(Indiv) %>%
  dplyr::mutate(POS = vcfR::getPOS(vcf.snp)) %>%
  dplyr::rename(acc = Indiv) %>%
  dplyr::left_join(., smpls, by = "acc") %>%
  dplyr::ungroup()


saveRDS(object = vcf.snp.mncnl.long,
        file = "data/derived_data/final_vcf_snps_long.RDS")



#........................................................
# Liftover to sequence call for fasta
#........................................................
mksngbp <- function(x){
  if (is.na(x)) {
    x <- "N"
  } else {
    x <- stringr::str_split_fixed(x, "/", n=2)[,1] # this [,1] works because I have forced to homozyg
  }
  return(x)
}

vcf.snp.mncnl.long$gt_GT_alleles_sbp <- sapply(vcf.snp.mncnl.long$gt_GT_alleles,
                                               mksngbp)

#........................................................
# write out final positions
#........................................................
# subset to seg sites now that we have coerced to no hets
pos <- unique( vcf.snp.mncnl.long$POS )
posgtcount <- vcf.snp.mncnl.long %>%
  dplyr::group_by(POS, gt_GT_alleles) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(!is.na(gt_GT_alleles))


# final seq sit positions
segsite <- unique(posgtcount$POS)
pos <- pos[pos %in% segsite]
saveRDS(object = pos, file = "data/derived_data/final_positions.rds")

#........................................................
# Make new fastas
#........................................................
seqlist <- vcf.snp.mncnl.long %>%
  dplyr::select(c("acc", "POS", "gt_GT_alleles_sbp"))
seqlist <- split(seqlist, factor(seqlist$acc))


mkmut <- function(x){
  mut <- pvp01
  mut[x$POS] <- x$gt_GT_alleles_sbp
  attr(mut, "name") <- unique(x$acc)
  return(mut)

}

finalfastas <- lapply(seqlist, mkmut)



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
dir.create(path = "data/fasta/", recursive = T)
seqinr::write.fasta(sequences = finalfastas,
                    names = names(finalfastas), file.out = "data/fasta/mtdna_anc.fa")


#......................
# ebro and drc
#......................
dir.create(path = "data/quick_look_fasta/", recursive = T)
drc_nha_ebro <- finalfastas[names(finalfastas) %in% c("D9U3K", "Ebro1944", "O6Y4K", "Q8J6O",
                                                      smpls$acc[smpls$vividregion == "NHA"])]
seqinr::write.fasta(sequences = drc_nha_ebro,
                    names = names(drc_nha_ebro), file.out = "data/quick_look_fasta/drc_nha_ebro.fa")

