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



#.................
# Class structure
#.................
plasm <- c(pv, po, pm, pk, pcyn, pchaub, pb)


#........................................................
# MSA
#........................................................






