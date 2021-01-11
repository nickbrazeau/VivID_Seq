#----------------------------------------------------
# Purpose of this script is to optimize and
# bootstrap phylogenetic trees
#----------------------------------------------------
# Internally, phangorn is removing duplicate sequences
# https://github.com/KlausVigo/phangorn/issues/16

library(ape)
library(phangorn)
#----------------------------------------------------
# Read In
#----------------------------------------------------
mtdna <- Biostrings::readDNAStringSet(filepath = "data/noclonefasta/unique_mtdna.fa")


#.................................
# read in data
#.................................
mtdna.ape <- ape::as.DNAbin(mtdna)
mtdna.dist.JC <- ape::dist.dna(mtdna.ape, model="JC69")


#.................................
# Phangorn
#.................................
mtdna.phangorn <- phangorn::as.phyDat(mtdna.ape)
tree.init <- ape::nj(mtdna.dist.JC)

#.................................
# JC Model
#.................................
JCfit.init <- phangorn::pml(tree.init, mtdna.phangorn)
fitJC <- phangorn::optim.pml(JCfit.init, model = "JC",
                             optNni=TRUE, optBf=TRUE, optQ=TRUE, optGamma=TRUE)


#.................................
# Bootstrap JC
#.................................
iters <- 1e3
JCbs <- phangorn::bootstrap.pml(fitJC, bs=iters,
                                 optNni=TRUE, optBf=TRUE, optQ=TRUE, optGamma=TRUE)


#.................................
# GTR Model
#.................................
fitGTR.init <- phangorn::pml(tree.init, mtdna.phangorn, k=4)
fitGTR <- phangorn::optim.pml(fitGTR.init, model="GTR",
                              optNni=TRUE, optBf=TRUE, optQ=TRUE, optGamma=TRUE)


#.................................
# Bootstrap GTR
#.................................
iters <- 1e3
GTRbs <- phangorn::bootstrap.pml(fitGTR, bs=iters,
                                 optNni=TRUE, optBf=TRUE, optQ=TRUE, optGamma=TRUE)


#.................................
# write out
#.................................
dir.create("all_smpl_tree_results")
saveRDS(object = fitJC, "all_smpl_tree_results/JCmodelfit.RDS")
saveRDS(object = JCbs, "all_smpl_tree_results/JC_bootstraptrees.RDS")

saveRDS(object = fitGTR, "all_smpl_tree_results/GTRfit.RDS")
saveRDS(object = GTRbs, "all_smpl_tree_results/GTR_bootstraptrees.RDS")
