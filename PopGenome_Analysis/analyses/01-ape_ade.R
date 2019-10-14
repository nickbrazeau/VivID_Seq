#----------------------------------------------------------------------------------------------------
# Purpose of this script is to perform analyses
# with the Ape, Adegenet, Ade4 suite
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(ape)
library(adegenet)

###############################################
####             Read in mtdt              ####
###############################################
smpls <- readxl::read_excel("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/vivid_seq_public_NGS.xlsx")

###############################################
####             Read in Seq               ####
###############################################
mtdna <- ape::read.FASTA(file = "data/derived_data/mtdna_anc.fa")

ama1haps <- Biostrings::DNAStringSet(ama1df$c_Consensus, use.names = T)
names(ama1haps) <- ama1df$h_popUID
ama1haps <- append(ama1haps, pfama1) # join ref
ama1haps <- as.matrix( as.DNAbin(ama1haps) ) # need matrix for boot phlyo, still dnabin










