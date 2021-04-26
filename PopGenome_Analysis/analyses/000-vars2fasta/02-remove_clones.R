#----------------------------------------------------
# Purpose of this script is to remove clones
# within VivID regions
#----------------------------------------------------
library(tidyverse)
library(Biostrings)
#----------------------------------------------------
# Read In
#----------------------------------------------------
mtdna <- Biostrings::readDNAStringSet(file = "data/fasta/mtdna_anc.fa")


smpls <- readxl::read_excel("../scrape_pubseqs/vivid_seq_public_NGS.xlsx") %>%
  dplyr::rename(smpls = acc) %>%
  dplyr::select(c("smpls", "country", "vividregion"))
fnlsmpls <- readRDS("data/derived_data/final_smpl_list.RDS")

fnlsmpls <- data.frame(smpls = fnlsmpls, stringsAsFactors = F)
fnlsmpls <- dplyr::left_join(fnlsmpls, smpls, by = "smpls")
fnlsmpls$country[is.na(fnlsmpls$country) & fnlsmpls$smpls %in%  c("D9U3K", "O6Y4K", "Q8J6O")] <- "CD"
fnlsmpls$vividregion[is.na(fnlsmpls$vividregion) & fnlsmpls$smpls %in%  c("D9U3K", "O6Y4K", "Q8J6O")] <- "AF"
# consider "apes" as own country
fnlsmpls$country[fnlsmpls$vividregion == "NHA"] <- "NHA"


# pull out lab samples which we will consider separately
labsmpls <- fnlsmpls[fnlsmpls$vividregion == "Lab",]
fnlsmpls <- fnlsmpls[fnlsmpls$vividregion != "Lab",]


#----------------------------------------------------
# Split List by Country so we can
# remove Clones by country
#----------------------------------------------------
fnlsmpls.split <- split(fnlsmpls$smpls, fnlsmpls$country)
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

mtdna.uniquehap.fulllist <- lapply(mtdna.list, findclones)


# save this out
saveRDS(object = mtdna.uniquehap.fulllist,
        file = "data/derived_data/mtdna_uniquehap_list.RDS")

# extract unique haps
mtdna.uniquehap.list <- purrr::map(mtdna.uniquehap.fulllist, "uniquehaps")

mtdna.unique <- NULL
for(i in 1:length(mtdna.uniquehap.list)){
  mtdna.unique <- append(mtdna.unique, mtdna.uniquehap.list[[i]])
}



#----------------------------------------------------
# Write out Clonal Corrections (sans lab) For Later Use
#----------------------------------------------------
dir.create("data/noclonefasta/")
Biostrings::writeXStringSet(x = mtdna.unique,
                            format = "fasta",
                            filepath = "data/noclonefasta/unique_mtdna.fa")


#----------------------------------------------------
# Write out Representative Sequences based on fastas
#----------------------------------------------------
# extract counts
mtdna_uniquehap_cnts <- purrr::map(mtdna.uniquehap.fulllist, "unihapcounts") %>%
  dplyr::bind_rows(., .id = "country")

# get representative
mtdna_rep_keep <- mtdna_uniquehap_cnts %>%
  dplyr::group_by(country) %>%
  dplyr::filter(hapuid_count == max(hapuid_count)) %>%
  dplyr::filter(!country %in% c("VN", "NHA")) %>% # see below
  dplyr::pull("smpls")

# random sample VN
set.seed(1234)
all_vn_smpls <- mtdna_uniquehap_cnts %>%
  dplyr::filter(country == "VN") %>%
  dplyr::pull("smpls")
VNsmpl <- sample(all_vn_smpls, 1)
# add back in NHA
NHAsmpls <- mtdna.uniquehap.fulllist$NHA$unihapcounts$smpls
# bring together
mtdna_rep_keep <- c(mtdna_rep_keep, VNsmpl, NHAsmpls)
# no make new small fasta
mtdna.rep <- mtdna[names(mtdna) %in% mtdna_rep_keep]

dir.create("data/rep_smpls_noclones/")
Biostrings::writeXStringSet(x = mtdna.rep,
                            format = "fasta",
                            filepath = "data/rep_smpls_noclones/rep_mtdna.fa")



#----------------------------------------------------
# Haplotype counts
#----------------------------------------------------
hapcount <- tibble::tibble(country = names(mtdna.list),
                           unihap = sapply(mtdna.uniquehap.list, length),
                           counthap = sapply(mtdna.list, length))
# write out
dir.create("tables")
readr::write_csv(hapcount, "tables/basic_hapcount.csv")
