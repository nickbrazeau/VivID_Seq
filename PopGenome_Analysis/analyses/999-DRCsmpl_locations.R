#..............................................................
# Purpose of this script is to plot where the three
# drc samples that we sequenced were from
#..............................................................
library(tidyverse)

#..............................................................
# DRC samples
#..............................................................
drcwgsmpls <- c("D9U3K", "O6Y4X", "Q8J6O")


#---------------------------------------------------------------------------------
# Pull in CD2013 Adult qPCR Results
#---------------------------------------------------------------------------------
pfpcr.adults <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Adults_Construction/Pf_alladults_v4.csv",
                                col_names = T) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::select(c("hivrecode_barcode", "pfldh", "fcq_mean")) %>%
  dplyr::rename(pfctmean = fcq_mean)

pvpcr.adults <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Adults_Construction/Pv_alladults_v2.csv",
                                col_names = T) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::select(c("hivrecode_barcode", "pv18s", "corrected_ct", "original_platemnum")) %>%
  dplyr::rename(pvctcrrct = corrected_ct)

pfpvpcr.adults <- dplyr::inner_join(pfpcr.adults, pvpcr.adults, by="hivrecode_barcode") %>%
  dplyr::rename(barcode = hivrecode_barcode,
                plate = original_platemnum) %>%
  dplyr::mutate(dataset = "adults")


#---------------------------------------------------------------------------------
# Pull in Kids qPCR Results
#---------------------------------------------------------------------------------
pfpvpcr.kids <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Children_Construction/PfPv_allchildren_v1.csv",
                                col_names = T) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(barcode = sh312,
                pvctcrrct = ct_corrected) %>%
  dplyr::mutate(dataset = "kids")


# note, these 7,250 kids are the under-5s versus the over-5s, match Mark's numbers perfectly
xtabs(~pfpvpcr.kids$pr_present, addNA = T)

# going to subset to just these under 5s and then
# drop to just important columns
pfpvpcr.kids <- pfpvpcr.kids %>%
  dplyr::filter(pr_present == 1) %>%
  dplyr::select("barcode", "pfldh", "pfctmean", "pv18s", "pvctcrrct",
                "plate", "dataset")



#.......................................
# Pull together kids and adults qpcr
#.......................................
qpcrdat <- dplyr::bind_rows(pfpvpcr.adults, pfpvpcr.kids)

#---------------------------------------------------------------------------------
# Join Kids and Adults qPCR data to DHS
#---------------------------------------------------------------------------------
# note, kids and adults have a duplicated barcode? i4b4r
# going to drop
qpcrdat <- qpcrdat %>%
  dplyr::filter(barcode != "i4b4r")
dt <- dplyr::inner_join(qpcrdat, arpr, by = "barcode")



#---------------------------------------------------------------------------------
# DHS spatial shapes/data
#---------------------------------------------------------------------------------

#spatial from the DHS -- these are cluster level vars
ge <- sf::st_as_sf(readRDS(file = "~/Documents/GitHub/VivID_KA_Compare/data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
colnames(ge) <- tolower(colnames(ge))
ge$adm1name <- gsub(" ", "-", ge$adm1name) # for some reason some of the char (like Kongo Central, lack a -), need this to match GADM
ge$adm1name <- gsub("Tanganyka", "Tanganyika", ge$adm1name) # DHS misspelled this province
# remove clusters that were missing from the DHS, see readme
ge <- ge %>%
  dplyr::rename(hv001 = dhsclust) # for easier merge with PR


dt <- left_join(dt, ge, by = "hv001")

# drop observations with missing geospatial data
dt.geo <- dt %>%
  dplyr::filter(latnum != 0 & longnum != 0) %>%
  dplyr::filter(!is.na(latnum) & !is.na(longnum))


