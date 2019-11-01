library(tidyverse)
library(RSelenium)


# Accessions
# For Cynomolgi Strain B acc - DRS000258 from PMC3759362
# For Cynomolgi Strain M acc - ERS001838 and ERS023609 from PMC5500898
# two above because it looks like they are technical replicates? Going to
# treat them as such
acc <- c("DRS000258", "ERS001838", "ERS023609")

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

# get rid of the non-PE illumina reads
enadf <- enadf %>%
  dplyr::filter(R2 != "")


readr::write_tsv(x=enadf,
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/Pcynomolgi/snakescrape_ENA/cynomolgi_scrape.tab.txt")
