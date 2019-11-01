#################################################################################
#### Purpose of this script is to find the ancestral alleles of the
#### P. vivax by determining which allele is ancestral (e.g. most conserved under the assumption of no recurrent mutation)
#### among the various Plasmodium species
#################################################################################

library(tidyverse)
library(Biostrings)
library(msa)

#........................................................
# Read in seqs
#........................................................
pv <- Biostrings::readDNAStringSet(filepath = "analyses/00-find_ancestral_alleles/sequences/PlasmoDB-45_PvivaxP01_Genome.fasta")
pv <- pv[16]
names(pv) <- "Pvivax_mtdna"

po <- readr::read_tsv(file = "analyses/00-find_ancestral_alleles/sequences/PocGH01_MIT_v1.embl",
                      col_names = F, skip = 4)
# for reasons I don't understand we concatenate a c on the beginning here, so need to drop it
po <- stringr::str_extract_all(po, "[actg]")[[1]]
po <- po[2:length(po)]
po <- Biostrings::DNAStringSet(paste(po, collapse = ""))
names(po) <- "Povale_mtdna"

pm <- Biostrings::readDNAStringSet(filepath = "analyses/00-find_ancestral_alleles/sequences/PlasmoDB-45_PmalariaeUG01_Genome.fasta")
pm <- pm[63]
names(pm) <- "Pmalariae_mtdna"

pk <- Biostrings::readDNAStringSet(filepath = "analyses/00-find_ancestral_alleles/sequences/PlasmoDB-45_PknowlesiH_Genome.fasta")
pk <- pk[164]
names(pk) <- "Pknowlesi_mtdna"

pcyn <- Biostrings::readDNAStringSet(filepath = "analyses/00-find_ancestral_alleles/sequences/PlasmoDB-45_PcynomolgiM_Genome.fasta")
pcyn <- pcyn[56]
names(pcyn) <- "Pcynomolgi_mtdna"

pchaub <- Biostrings::readDNAStringSet(filepath = "analyses/00-find_ancestral_alleles/sequences/PlasmoDB-45_Pchabaudichabaudi_Genome.fasta")
pchaub <- pchaub[16]
names(pchaub) <- "Pchabaudi_mtdna"

pb <- Biostrings::readDNAStringSet(filepath = "analyses/00-find_ancestral_alleles/sequences/PlasmoDB-45_PbergheiANKA_Genome.fasta")
pb <- pb[21]
names(pb) <- "Pberghei_mtdna"


pf <- Biostrings::readDNAStringSet(filepath = "analyses/00-find_ancestral_alleles/sequences/PlasmoDB-45_Pfalciparum3D7_Genome.fasta")
pf <- pf[16]
names(pf) <- "Pfalciparum_mtdna"


#.................
# Class structure
#.................
plasm <- c(pv, po, pm, pk, pcyn, pchaub, pb, pf)
Biostrings::writeXStringSet(x=plasm,
                            filepath = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/PopGenome_Analysis/analyses/00-find_ancestral_alleles/Plasmodium_MtDNA.fasta",
                            format = "fasta")

#........................................................
# MARS and MSA
#........................................................
# system("bash runMars.sh")
# both used default settings
# Muscle version 3.8.31 from plasm.mars.msa@version
plasm.mars <- Biostrings::readDNAStringSet(filepath = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/PopGenome_Analysis/analyses/00-find_ancestral_alleles/Plasmodium_MtDNA.mars.fasta")
plasm.mars.msa <- msa::msaMuscle(plasm.mars)

# make sure Pv doesn't get rotated
plasm.mars[1][[1]] == plasm[1][[1]]

sum( seqinr::s2c(as.character( plasm.mars[1][[1]] )) !=
       seqinr::s2c(as.character( plasm[1][[1]] ))
)

#........................................................
# Get Ancestral Allele
#........................................................
ancestralfa <- rep(NA, unique(plasm.mars.msa@unmasked@ranges@width))
msamat <- as.matrix( plasm.mars.msa@unmasked )
# reorder
pv <- which(rownames( msamat ) == "Pvivax_mtdna")
seq <- 1:nrow(msamat)
seq <- seq[! 1:nrow(msamat) %in% pv ]

msamat.ord <- msamat[c(pv,seq), ]

# start walking along genome
for(i in 1:ncol(msamat.ord)){
  if(msamat.ord[1,i] == "-"){ # if pv has a gap, keep gap
    ancestralfa[i] <- "-"
  } else if(all(msamat.ord[2:nrow(msamat.ord),i] == "-")) {
    ancestralfa[i] <- "N"
  }
  else{
    alleletab <- table(msamat.ord[,i])
    allele <- names(alleletab)[alleletab == max(alleletab)]
    if(length(allele) > 1 | any(allele == "-")){ # not enough support to call ancestral state
      ancestralfa[i] <- "N"
    } else{
      ancestralfa[i] <- allele
    }
  }
}

ancestralfa <- ancestralfa[ancestralfa != "-"]
ancestralfa <- Biostrings::DNAStringSet( paste(ancestralfa, collapse = "") )
names(ancestralfa) <- "ancestral_allele"

#........................................................
# Sanity Check
#........................................................
orig <- plasm[1]
names(orig) <- "orig"
mars <- plasm.mars[1]
names(mars) <- "mars"
test <- ancestralfa
names(ancestralfa) <- "ancestral"
out <- c(orig, mars, test)
writeXStringSet(out, filepath="~/Desktop/sanitycheck.fasta")

# looks reasonsable
writeXStringSet(ancestralfa, filepath="~/Documents/GitHub/VivID_Seq/PopGenome_Analysis/analyses/00-find_ancestral_alleles/ancestral.fa")










