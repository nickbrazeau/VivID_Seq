library(tidyverse)
library(genoPlotR)
library(Biostrings)
library(ape)

#...................
# Read in seqs
#...................
pv <- Biostrings::readDNAStringSet(filepath = "~/Documents/MountPoints/mountIDEEL/resources/genomes/Pvivax/genomes/PvP01.fasta")
pv <- pv[grepl("MIT", pv@ranges)]
names(pv) <- "Pvivax_mtdna"
po <- Biostrings::readDNAStringSet(filepath = "~/Documents/GitHub/VivID_Seq/find_ancestral_alleles/sequences/Povale_contigs.fasta.gz")
po <- po[grepl("MIT", po@ranges)]
names(po) <- "Povale_mtdna"
pm <- Biostrings::readDNAStringSet(filepath = "~/Documents/GitHub/VivID_Seq/find_ancestral_alleles/sequences/Pmalariae_contigs.fasta.gz")
pm <- pm[grepl("MIT", pm@ranges)]
names(pm) <- "Pmalariae_mtdna"

pk <- readr::read_tsv("~/Documents/GitHub/VivID_Seq/find_ancestral_alleles/sequences/PKNH_MIT_v2.embl",
                      comment = c("FT"), skip = 3, col_names = T) 
pk <- stringr::str_extract_all(pk, "[actg]", simplify = F)[[1]]
# for reasons I don't understand there is a c being appeneded here
pk <- pk[2:length(pk)]
pk <- Biostrings::DNAStringSet(paste(pk, collapse = ""))
names(pk) <- "Pknowlesi_mtdna"


pc <- readr::read_tsv("~/Documents/GitHub/VivID_Seq/find_ancestral_alleles/sequences/PcyM_MT.embl",
                      comment = c("FT"), skip = 5, col_names = T) 
pc <- stringr::str_extract_all(pc, "[actg]", simplify = F)[[1]]
# for reasons I don't understand there is a c being appeneded here
pc <- pc[2:length(pc)]
pc <- Biostrings::DNAStringSet(paste(pc, collapse = ""))
names(pc) <- "Pcynomolgi_mtdna"

#...................
# MSA
#...................








temp = readr::read_tsv("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/ENA_master_acc_download_map_PE.tab.txt")
long <- temp %>% 
  dplyr::select(-c("acc")) %>% 
  unlist(.)

readr::write_tsv(x=data.frame(long), path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/ENA_master_acc_download_map_PE-LONG.tab.txt",
                 col_names = F)
