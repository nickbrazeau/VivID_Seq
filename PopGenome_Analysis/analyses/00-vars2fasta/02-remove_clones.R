#----------------------------------------------------
# Purpose of this script is to remove clones
# within VivID regions
#----------------------------------------------------
library(tidyverse)
library(Biostrings)
#----------------------------------------------------
# Read In
#----------------------------------------------------
mtdna <- Biostrings::readDNAStringSet(file = "~/Documents/GitHub/VivID_Seq/PopGenome_Analysis/data/fasta/mtdna_anc.fa")


smpls <- readxl::read_excel("~/Documents/GitHub/VivID_Seq/scrape_pubseqs/vivid_seq_public_NGS.xlsx") %>%
  dplyr::rename(smpls = acc) %>%
  dplyr::select(c("smpls", "country", "vividregion"))
fnlsmpls <- readRDS("~/Documents/GitHub/VivID_Seq/PopGenome_Analysis/data/derived_data/final_smpl_list.RDS")

fnlsmpls <- data.frame(smpls = fnlsmpls, stringsAsFactors = F)
fnlsmpls <- dplyr::left_join(fnlsmpls, smpls, by = "smpls")
fnlsmpls$country[is.na(fnlsmpls$country) & fnlsmpls$smpls %in%  c("D9U3K", "O6Y4K", "Q8J6O")] <- "CD"
fnlsmpls$vividregion[is.na(fnlsmpls$vividregion) & fnlsmpls$smpls %in%  c("D9U3K", "O6Y4K", "Q8J6O")] <- "AF"

#----------------------------------------------------
# Split List by Country so we can
# remove Clones by country
#----------------------------------------------------
fnlsmpls.split <- split(fnlsmpls$smpls, fnlsmpls$country)
fnlsmpls.split <- append(fnlsmpls.split, list(anc = "Pcynomolgi"))
mtdna.list <- lapply(fnlsmpls.split, function(x) return( mtdna[ names(mtdna) %in% x ] ))

#----------------------------------------------------
# Find Clones
# (not considering missing data)
#----------------------------------------------------
findclones <- function(countrymtdna){

  unihaps <- unique(countrymtdna)
  unihapcounts <- lapply(unihaps, function(x){
    ret <- length(which(countrymtdna == x))
    return(ret)
  })

  unihapcounts <- tibble::tibble(smpls = names(unihapcounts),
                                 hapuid_count = unlist(unihapcounts))

  ret <- list(uniquehaps = unihaps,
              unihapcounts = unihapcounts)
  return(ret)

}

mtdna.uniquehap.list <- lapply(mtdna.list, findclones)



#----------------------------------------------------
# Write out Clonal Corrections For Later Use
#----------------------------------------------------
saveRDS(object = mtdna.uniquehap.list,
        file = "data/derived_data/mtdna_uniquehap_list.RDS")

mtdna.uniquehap.list <- purrr::map(mtdna.uniquehap.list, "uniquehaps")

mtdna.unique <- NULL
for(i in 1:length(mtdna.uniquehap.list)){
  mtdna.unique <- append(mtdna.unique, mtdna.uniquehap.list[[i]])
}

dir.create("data/noclonefasta/")
Biostrings::writeXStringSet(x = mtdna.unique,
                            format = "fasta",
                            filepath = "data/noclonefasta/unique_mtdna.fa")
