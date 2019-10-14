library(tidyverse)
library(rentrez)
setwd("/proj/ideel/meshnick/users/NickB/Projects/VivID_Seq/scrape_pubseqs")
ncbiapikey <- readr::read_tsv("ncbi_api_key_forsanger.txt", col_names = F)
rentrez::set_entrez_key(unname(unlist(ncbiapikey)))

mtsanger <- readxl::read_excel("vivid_seq_public_sangerseq.xlsx")


## Query
seqacc <- tibble::tibble(acc = mtsanger$acc)
search <- lapply(1:nrow(seqacc), function(x) return(list())) # nested list
for(query in 1:nrow(seqacc)){

  acc.search <- rentrez::entrez_search(db= "nucleotide",
                                       seqacc$acc[query],
                                       retmax = 5,
                                       use_history = T)
  search[[query]] <- acc.search

}

## Get the IDs that we just found
seqacc$id <- lapply(search, function(x){x$ids})

## Scrape/Fetch
# https://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
fasta_fetch <- function(id){
  ret <- rentrez::entrez_fetch(db="nucleotide", rettype="fasta", id = id)
  Sys.sleep(60)
  return(ret)
}


seqacc$fasta <- NA
for(i in 1:nrow(seqacc)){
  seqacc$fasta[i] <- fasta_fetch(seqacc$id[i])
}

dir.create( "/pine/scr/n/f/nfb/Projects/VivID_Seq/public_sanger_seqs/", recursive = T)
saveRDS(seqacc, "/pine/scr/n/f/nfb/Projects/VivID_Seq/public_sanger_seqs/global_sangerseq_seqacc.rds")
seqinr::write.fasta(sequences = outfastas,
                    names = seqacc$acc,
                    file.out = "/pine/scr/n/f/nfb/Projects/VivID_Seq/public_sanger_seqs/global_sangerseq_vivid.fasta")





