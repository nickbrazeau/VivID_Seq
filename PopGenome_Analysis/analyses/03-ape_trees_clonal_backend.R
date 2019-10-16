library(ape)
#----------------------------------------------------
# Read In
#----------------------------------------------------
mtdna <- Biostrings::readDNAStringSet(filepath = "data/noclonefasta/unique_mtdna.fa")
# since we are not estimating the molecular clock
mtdna <- mtdna[!names(mtdna) %in% "Anc"]


#.................................
# read in data
#.................................
mtdna.unique.ape <- ape::as.DNAbin(mtdna)
mtdna.dist.JC <- ape::dist.dna(mtdna.unique.ape, model="JC69")


#.................................
# Phangorn
#.................................
mtdna.unique.phangorn <- phangorn::as.phyDat(mtdna.unique.ape)
tree.init <- nj(mtdna.dist.JC)

#.................................
# JC Model
#.................................
JCfit.init <- phangorn::pml(tree.init, mtdna.unique.phangorn)
fitJC <- optim.pml(JCfit.init, model = "JC",
                   optNni=TRUE, optBf=TRUE, optQ=TRUE)


#.................................
# Bootstrap JC
#.................................
iters <- 1e3
JCbs <- phangorn::bootstrap.pml(fitJC, bs=iters,
                                 optNni=TRUE, optBf=TRUE, optQ=TRUE, optGamma=TRUE)


#.................................
# GTR Model
#.................................
fitGTR.init <- phangorn::pml(tree.init, mtdna.unique.phangorn, k=4)
fitGTR <- optim.pml(fitGTR.init, model="GTR",
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
dir.create("clonal_smpl_tree_results")
saveRDS(object = fitJC, "clonal_smpl_tree_results/JCmodelfit.RDS")
saveRDS(object = JCbs, "clonal_smpl_tree_results/JC_bootstraptrees.RDS")

saveRDS(object = fitGTR, "clonal_smpl_tree_results/GTRfit.RDS")
saveRDS(object = GTRbs, "clonal_smpl_tree_results/GTR_bootstraptrees.RDS")
