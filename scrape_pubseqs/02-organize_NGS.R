library(tidyverse)


#...........................................................................................
# Organize PE for Alignment
#............................................................................................
pe <- readr::read_tsv("~/Documents/MountPoints/mountMeshnick/Projects/VivID_Seq/scrape_pubseqs/ENA_master_acc_download_map_PE.tab.txt")
pe.long <- pe %>%
  dplyr::mutate(SRR1 = basename(R1),
                SRR2 = basename(R2)) %>%
  dplyr::select(-c("R1", "R2")) %>%
  tidyr::gather(., key = "srr", value = "fastq", 2:3) %>%
  dplyr::select(-c("srr"))

pe.paths <- tibble::tibble(path = list.files(path = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/public_pe_seqs/",
                                              pattern = "fastq.gz",
                                              full.names = T)) %>%
  dplyr::mutate(fastq = basename(path))

pe.long.paths <- dplyr::left_join(pe.long, pe.paths, by = "fastq")



#...............................
# Make Symlink Architecture
#................................
symlink_architecture <- pe.long.paths %>%
  magrittr::set_colnames(c("smpl", "fastq", "from")) %>%
  dplyr::mutate(to = paste0(smpl, "/", fastq)) %>%
  dplyr::select(-c("fastq")) %>%
  dplyr::mutate(from = gsub("/Users/nbrazeau/Documents/MountPoints/mountedScratchLL/",
                            "/pine/scr/n/f/nfb/",
                            from))


readr::write_tsv(x = symlink_architecture,
                 path = "~/Documents/MountPoints/mountMeshnick/Projects/VivID_Seq/wgs_pe_improved_global/symlink_architecture.tab.txt",
                 col_names = F)

#...............................
# Make Run Map
#................................
globalvivid_run_map <- pe.long.paths %>%
  dplyr::select(-c("path")) %>%
  dplyr::mutate(fastq = stringr::str_split_fixed(fastq, "_", n=2)[,1]) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::mutate(x = ".")

readr::write_tsv(x = globalvivid_run_map,
                 path = "~/Documents/MountPoints/mountMeshnick/Projects/VivID_Seq/wgs_pe_improved_global/globalvivid_run_map.tab.txt",
                 col_names = F)




#...........................................................................................
# Organize SE for Alignment
#............................................................................................
se <- readr::read_tsv("~/Documents/MountPoints/mountMeshnick/Projects/VivID_Seq/scrape_pubseqs/ENA_master_acc_download_map_SE.tab.txt")

se.long <- se %>%
  dplyr::mutate(SRR1 = basename(R1)) %>%
  dplyr::select(-c("R1")) %>%
  tidyr::gather(., key = "srr", value = "fastq", 2) %>%
  dplyr::select(-c("srr"))

se.paths <- tibble::tibble(path = list.files(path = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/public_se_seqs/",
                                             pattern = "fastq.gz",
                                             full.names = T)) %>%
  dplyr::mutate(fastq = basename(path))

# Note, SYpte56 (SRS3371818 - has pacbio reads SRR7255036 and illumina reads SRR7255037)
# Note, SYptt43 (SRS3371817 - has pacbio reads SRR7255038 and illumina reads SRR7255039)
# have to drop the pacbio reads
se.paths <- se.paths[ !grepl("SRR7255036|SRR7255038", se.paths$fastq), ]

se.long.paths <- dplyr::left_join(se.long, se.paths, by = "fastq") %>%
  dplyr::filter(!is.na(path))

# Note, India samples SRS805922, SRS807702, SRS807711, SRS807712, SRS805942, SRS805943, SRS807544, SRS807701
# have both PE and SE reads. Drop all SE reads (as we will not combine these) due to different error mode

se.long.paths <- se.long.paths[!se.long.paths$acc %in% pe.long.paths$acc, ]

#......................
# manually add in Ebro from study authors
#......................
se.long.paths <- tibble::tibble(acc = "Ebro1944",
                                fastq = "Ebro1944.fastq.gz",
                                path = "/pine/scr/n/f/nfb/Projects/VivID_Seq/public_se_seqs/Ebro1944.fastq.gz")

#...............................
# Make Symlink Architecture (overkill)
#................................
symlink_architecture <- se.long.paths %>%
  magrittr::set_colnames(c("smpl", "fastq", "from")) %>%
  dplyr::mutate(to = paste0(smpl, "/", fastq)) %>%
  dplyr::select(-c("fastq")) %>%
  dplyr::mutate(from = gsub("/Users/nickbrazeau/Documents/MountPoints/mountedScratchLL/",
                            "/pine/scr/n/f/nfb/",
                            from))

readr::write_tsv(x = symlink_architecture,
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/wgs_se_improved_global/symlink_architecture.tab.txt",
                 col_names = F)

#...............................
# Make Run Map
#................................
globalvivid_run_map <- se.long.paths %>%
  dplyr::select(-c("path")) %>%
  dplyr::mutate(fastq = gsub(".fastq.gz", "", fastq)) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::mutate(x = ".")

readr::write_tsv(x = globalvivid_run_map,
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/wgs_se_improved_global/globalvivid_run_map.tab.txt",
                 col_names = F)






