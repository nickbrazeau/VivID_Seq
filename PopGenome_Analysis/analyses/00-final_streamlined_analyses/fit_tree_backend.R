## .................................................................................
## Purpose: From previous results, I know GTR model is best. Going to
## run out on cluster for time so I can do more boostrap
##
## .................................................................................
set.seed(48)
library(tidyverse)
library(phangorn)
library(ape)
library(Biostrings)
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


#............................................................
# bootstrap tree
#...........................................................
mtdna.ape <- ape::as.DNAbin(mtdna_uniq)
tree.init <- ape::nj(ape::dist.dna(mtdna.ape, model="N"))

# Phangorn
mtdna.phangorn <- phangorn::as.phyDat(mtdna.ape)
# GTR Model
fitGTR.init <- phangorn::pml(tree.init, mtdna.phangorn, k=4)
fitGTR <- phangorn::optim.pml(fitGTR.init, model = "GTR",
                              rearrangement = "stochastic",
                              optNni=TRUE, optBf=TRUE, optQ=TRUE, optGamma=TRUE)

# Bootstrap GTR which is better model
iters <- 1e3
GTRbs <- phangorn::bootstrap.pml(fitGTR, bs=iters,
                                 optNni = TRUE,
                                 rearrangement = "stochastic")

dir.create("analyses/00-final_streamlined_analyses/tree_results_sm", recursive = T)
saveRDS(fitGTR, "analyses/00-final_streamlined_analyses/tree_results_sm/tree_ml.rds")
saveRDS(GTRbs, "analyses/00-final_streamlined_analyses/tree_results_sm/tree_boot.rds")
