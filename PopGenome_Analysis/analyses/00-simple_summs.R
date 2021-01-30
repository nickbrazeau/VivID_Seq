## .................................................................................
## Purpose: Simple Script to bring together final results
##          to include for phylogenetic analysis
##
## Notes: This is the only file that matters in this analysis repo.
##        The other RMD scripts are depreciated but kept along for
##        those interested in the additional data
## .................................................................................
library(tidyverse)
library(PopGenome)
library(Biostrings)

#............................................................
###### Section 1: Quick stats #####
#...........................................................
#......................
# in data
#......................
mtdt <- readRDS("data/derived_data/mtdt.RDS")
# long vcf
vcf.snp.ss.long_gt <- readRDS("data/derived_data/final_vcf_snps_long.RDS")
vcf.snp.ss.long_gt <- vcf.snp.ss.long_gt %>%
  dplyr::left_join(., y = mtdt, by = c("iid", "country", "vividregion")) %>%
  dplyr::filter(passed == "Y") %>%
  dplyr::mutate(POS = factor(POS, levels = c(1:5989), ordered = T))

# all dna
mtdna_all <- Biostrings::readDNAStringSet(filepath = "data/fasta/mtdna_anc.fa")
# no clones
mtdna_uniq <- Biostrings::readDNAStringSet(filepath = "data/noclonefasta/unique_mtdna.fa")


#......................
# basic popgen statistics
#......................
mtdna_uniq
mtdna_all
# get seg sites
seqmat <-  tibble::as_tibble( t(as.matrix(mtdna_all) ) )
seqmat <- cbind.data.frame(POS = 1:nrow(seqmat), seqmat)
segsites <- apply(seqmat[,2:ncol(seqmat)], 1, function(x){
  ret <- unique(x)
  ret <- ret[ret != "N"]
  return(ret)
})
sum(sapply(segsites, length) != 1)

# define countries
mtdt <- mtdt %>%
  dplyr::mutate(effcountry = ifelse(vividregion == "Lab", "Lab",
                                    ifelse(vividregion == "NHA", "NHA",
                                           country)))
# cynomolgi
mtdt <- mtdt %>%
  dplyr::mutate(effcountry = ifelse(smpls == "ERS001838", "PcyM", effcountry))


# read in and set up popgenome object
fa_popgen <- PopGenome::readData(path = "data/fasta/",
                                 format = "fasta",
                                 include.unknown = T,
                                 progress_bar_switch=T,
                                 FAST=F,
                                 big.data=F)
rgn <- split(mtdt$smpls, factor(mtdt$effcountry))
fa_popgen <- PopGenome::set.populations(fa_popgen,
                                        rgn)
# set Cynmologi as outgroup
fa_popgen <- PopGenome::set.outgroup(fa_popgen,
                                     c("ERS001838")) #PcyM-strain
get.sum.data(fa_popgen)

#......................
# w/in diversity
#......................
# must be calculated first before diversity and neutrality
fa_popgen <- PopGenome::F_ST.stats(fa_popgen)
fa_popgen <- PopGenome::diversity.stats.between(fa_popgen)
fa_popgen <- PopGenome::neutrality.stats(fa_popgen)
ret <- PopGenome::get.diversity(fa_popgen, between = F)
popn <- gsub(" ", "", rownames(ret))
tbout <- lapply(ret, as.data.frame) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(country = names(fa_popgen@populations)) %>%
  dplyr::select(c("country", "nuc.diversity.within")) %>%
  dplyr::mutate(nuc.diversity.within = nuc.diversity.within/fa_popgen@n.sites) # PER DOC -- The nucleotide diversity (average pairwise nucleotide differences) within and between the populations. Have to be divided by the slot GENOME.class@n.sites to obtain diversity per site (see also diversity.stats).

# write out
dir.create("tables")
readr::write_csv(tbout, path = "tables/basic_popgen_summstat.csv")

#............................................................
# basic hammings distance
#...........................................................
hammdist <- ape::dist.dna(ape::as.DNAbin(mtdna_all), model = "N")
zerodist <- hammdist %>%
  broom::tidy(.) %>%
  dplyr::filter(distance == 0) %>%
  dplyr::filter(item1 %in% c("D9U3K", "O6Y4K", "Q8J6O") |
                  item2 %in% c("D9U3K", "O6Y4K", "Q8J6O"))

zerodist <- unique(c(as.character(zerodist$item1), as.character(zerodist$item2)))
zerodist <- zerodist[!zerodist %in% c("D9U3K", "O6Y4K", "Q8J6O")]
# get smpls
zerodist_smpls <- mtdt %>%
  dplyr::filter(smpls %in% zerodist)
# look to see if all same haplotype for TH
mtdna.uniquehap.list <- readRDS("data/derived_data/mtdna_uniquehap_list.RDS")
mtdna.uniquehap.list[["TH"]]$unihapcounts[
  mtdna.uniquehap.list[["TH"]]$unihapcounts$smpls %in% zerodist_smpls$smpls,]

#......................
# nha and ebro dist
#......................
nhasmpls <- mtdt %>%
  dplyr::filter(vividregion == "NHA") %>%
  dplyr::pull("smpls")
hammdist %>%
  broom::tidy(.) %>%
  dplyr::filter(item1 %in% c("D9U3K", "O6Y4K", "Q8J6O") |
                  item2 %in% c("D9U3K", "O6Y4K", "Q8J6O")) %>%
  dplyr::filter(item1 %in% c(nhasmpls, "Ebro1944") |
                  item2 %in% c(nhasmpls, "Ebro1944"))

#......................
# all
#......................
alldrcpairs <- hammdist %>%
  broom::tidy(.) %>%
  dplyr::filter(item1 %in% c("D9U3K", "O6Y4K", "Q8J6O") |
                  item2 %in% c("D9U3K", "O6Y4K", "Q8J6O")) %>%
  dplyr::filter(!(item1 == "ERS001838" |
                    item2 == "ERS001838")) # drop PyCm
summary(alldrcpairs)

#............................................................
###### Section 2: GT Viz #####
#...........................................................
# find seg sites from fasta
seqmat <-  tibble::as_tibble( t(as.matrix(mtdna_all) ) )
seqmat <- cbind.data.frame(POS = 1:nrow(seqmat), seqmat)
segsites <- apply(seqmat[,2:ncol(seqmat)], 1, function(x){ return( length(unique(x)) > 1 ) })

# take to long
seqmat.long <- seqmat %>%
  dplyr::filter(segsites) %>%
  tidyr::gather(., key = "smpls", value = "GT", 2:ncol(.)) %>%
  dplyr::mutate(GT = ifelse(GT == "N", NA, GT)) %>%
  dplyr::mutate(GT = factor(GT)) %>%
  dplyr::left_join(., mtdt, by = "smpls") %>%
  dplyr::mutate(POS = factor(POS, levels = c(1:5989), ordered = T))


# make plot
gtplot <- seqmat.long %>%
  ggplot() +
  geom_tile(aes(x = smpls, y = POS, fill = GT)) +
  scale_fill_manual("Alleles", values = c("#4daf4a", "#377eb8", "#ffff33", "#e41a1c")) +
  ylab("POS") +
  facet_wrap(~effcountry, scales = "free_x", nrow = 1, strip.position = "top") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_blank())

jpeg("figures/consensus_haplotypes_vivax.jpg",
     width = 11, height = 8, unit = "in", res = 500)
plot(gtplot)
graphics.off()
