ncbiapikey <- readr::read_tsv("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/ncbi_api_key_forsanger.txt", col_names = F)
set_entrez_key(unname(unlist(ncbiapikey)))
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
# https://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
fasta_fetch <- function(id){
  out <- tryCatch(
    {
      rentrez::entrez_fetch(db="nucleotide", rettype="fasta", id = id)
      Sys.sleep(10)
    },
    error=function(cond) {
      message(paste("ID does not seem to exist:", id))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("ID caused a warning:", id))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      message(paste("Processed id:", id))
    }
  )    
  return(out)
}
  



seqacc$fasta <- purrr::map(seqacc$id, fasta_fetch)



outfastas <- purrr::map(seqacc, `[[`, "fasta")
dir.create("~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/vividseq_sanger_publicseqs/",
           recursive = T)
write.table(x=outfastas, file="~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/vividseq_sanger_publicseqs/vividseq_sanger_publicseqs.fa", 
            quote = F, col.names = F, row.names = F, sep = "")





