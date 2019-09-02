library(tidyverse)


#...........................................................................................
# Organize PE for Alignment
#............................................................................................
pe <- readr::read_tsv("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/ENA_master_acc_download_map_PE.tab.txt")
pe.long <- pe %>% 
  dplyr::mutate(SRR1 = basename(R1),
                SRR2 = basename(R2)) %>% 
  dplyr::select(-c("R1", "R2")) %>% 
  tidyr::gather(., key = "srr", value = "fastq") %>% 
  dplyr::select(-c("srr"))

pe.paths <- tibble::tibble(path = list.files(path = "",
                                              pattern = "fastq.gz",
                                              full.names = T)) %>% 
  dplyr::mutate(fastq = basename(path))

#...............................
# Make Symlink Architecture 
#................................


#...............................
# Make Run Map
#................................

pe.run.map <- dplyr::left_join(pe.long, pe.paths) %>% 
  dplyr::mutate(read = ifelse(stringr::str_detect(fastq, "_1.fastq.gz"), "R1", "R2")
                ) %>% 
  tidyr::spread(., key = "read", value = "path")




#...........................................................................................
# Organize SE for Alignment
#............................................................................................
se <- readr::read_tsv("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/ENA_master_acc_download_map_SE.tab.txt")
