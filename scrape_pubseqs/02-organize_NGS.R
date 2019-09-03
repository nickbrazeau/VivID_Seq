library(tidyverse)


#...........................................................................................
# Organize PE for Alignment
#............................................................................................
pe <- readr::read_tsv("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/ENA_master_acc_download_map_PE.tab.txt")
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
  dplyr::mutate(from = gsub("/Users/nickbrazeau/Documents/MountPoints/mountedScratchLL/", 
                            "/proj/ideel/meshnick/users/NickB/", 
                            from))
    
readr::write_tsv(x = symlink_architecture, 
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/wgs_pe_improved_global/symlink_architecture.tab.txt",
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
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/wgs_pe_improved_global/globalvivid_run_map.tab.txt",
                 col_names = F)  






#...........................................................................................
# Organize SE for Alignment
#............................................................................................
se <- readr::read_tsv("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/ENA_master_acc_download_map_SE.tab.txt")
