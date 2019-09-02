library(tidyverse)
library(RSelenium)


# read in metadata
smpls <- readxl::read_excel("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/vivid_seq_public_NGS.xlsx")
acc <- smpls$acc
# drive 
remDr <- RSelenium::rsDriver(verbose = F, port = 3644L, browser = "chrome")  # instantiate remote driver to connect to Selenium Server

for(smpl in 1:length(acc)){
  # assign the client to a new variable
  RD <- remDr$client
  RD$navigate(paste0("https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=", 
                     acc[smpl],
                     "&result=read_run&fields=fastq_ftp&download=txt")
  )

  Sys.sleep(5) 

}

RD$close()
remDr$closeServer

enafiles <- list.files(path = "~/Downloads/", pattern = ".txt", full.names = T)
readena <- function(path){
  ret <- readr::read_tsv(path) %>% 
    dplyr::mutate(acc = gsub(".txt", "", basename(path)), 
                  R1 = stringr::str_split_fixed(fastq_ftp, ";", n=2)[,1],
                  R2 = stringr::str_split_fixed(fastq_ftp, ";", n=2)[,2]) %>% 
    dplyr::select(c("acc", "R1", "R2"))
  
  return(ret)
  
}

enadf <- lapply(enafiles, readena) %>% 
  dplyr::bind_rows()

#...............
# Single End
#...............
enadf.se <- enadf %>% 
  dplyr::filter(R2 == "") %>% 
  dplyr::select(-c("R2"))

#...............
# Paired End
#...............
enadf.pe <- enadf %>% 
  dplyr::filter(R2 != "") 

# note PMC6032790 has an extra fastq (extra file) that isn't PE ? Going to subset just to R1/R2
enadf.pe.norm <- enadf.pe %>% 
  dplyr::filter(!grepl(";", R2))

enadf.pe.extra <- enadf.pe %>% 
  dplyr::filter(grepl(";", R2)) %>% 
  dplyr::mutate(temp = R2,
                R1 = NA, 
                R2 = NA,
                R1 = stringr::str_split_fixed(temp, ";", n=2)[,1],
                R2 = stringr::str_split_fixed(temp, ";", n=2)[,2]
                ) %>% 
  dplyr::select(-c("temp"))

enadf.pe.cleaned <- rbind.data.frame(enadf.pe.norm, enadf.pe.extra) 

readr::write_tsv(x=enadf.pe.cleaned, 
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/ENA_master_acc_download_map_PE.tab.txt")

readr::write_tsv(x=enadf.se, 
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/ENA_master_acc_download_map_SE.tab.txt")

