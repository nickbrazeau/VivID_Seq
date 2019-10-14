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
  dplyr::select(c("smpls", "vividregion"))
fnlsmpls <- readRDS("~/Documents/GitHub/VivID_Seq/PopGenome_Analysis/data/derived_data/final_smpl_list.RDS")

fnlsmpls <- data.frame(smpls = fnlsmpls, stringsAsFactors = F)
fnlsmpls <- dplyr::left_join(fnlsmpls, smpls, by = "smpls")
fnlsmpls$vividregion[is.na(fnlsmpls$vividregion)] <- "DRC"

# split final samples
fnlsmpls.split <- split(fnlsmpls$smpls, fnlsmpls$vividregion)
fnlsmpls.split <- append(fnlsmpls.split, list(anc = "Anc"))

mtdna.list <- lapply(fnlsmpls.split, function(x) return( mtdna[ names(mtdna) %in% x ] ))

# fine clones (not considering missing data)
mtdna.uniquehap.list <- lapply(mtdna.list, unique)


#----------------------------------------------------
# Write out Clonal Corrections For Later Use
#----------------------------------------------------
saveRDS(object = mtdna.uniquehap.list,
        file = "data/derived_data/mtdna_uniquehap_list.RDS")

mtdna.unique <- NULL
for(i in 1:length(mtdna.uniquehap.list)){
  mtdna.unique <- append(mtdna.unique, mtdna.uniquehap.list[[i]])
}

dir.create("data/noclonefasta/")
Biostrings::writeXStringSet(x = mtdna.unique,
                            format = "fasta",
                            filepath = "data/noclonefasta/unique_mtdna.fa")
