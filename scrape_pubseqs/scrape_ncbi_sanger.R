library(tidyverse) 
library(rentrez)

mtsanger <- readxl::read_excel("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/vivid_seq_public_sangerseq.xlsx")


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
fasta_fetch <- function(id){
  f <- rentrez::entrez_fetch(db="nucleotide", rettype="fasta", id = id)
  Sys.sleep(10)
  return(f)
}

seqacc$fasta <- NA
for(i in 1:length(seqacc$id)){
  seqacc$fasta[i] <- fasta_fetch(seqacc$id[i])
} 


outfastas <- purrr::map(seqacc, `[[`, "fasta")
dir.create("~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/vividseq_sanger_publicseqs/",
           recursive = T)
write.table(x=outfastas, file="~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/vividseq_sanger_publicseqs/vividseq_sanger_publicseqs.fa", 
            quote = F, col.names = F, row.names = F, sep = "")





