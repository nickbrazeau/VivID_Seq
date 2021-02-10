####################################################################################
## Purpose: Haplotype Network Based
##
## Notes:
####################################################################################
set.seed(48)
library(ape)
library(tidyverse)
#............................................................
# read in
#...........................................................
mtdna_uniq <- Biostrings::readDNAStringSet(filepath = "data/noclonefasta/unique_mtdna.fa")
# mtdt
mtdt <- readRDS("data/derived_data/mtdt.RDS") %>%
  dplyr::filter(smpls %in% names(mtdna_uniq)) %>%
  dplyr::mutate(country = ifelse(vividregion == "NHA", "NHA",
                                 ifelse(iid == "Ebro1944", "Ebro1944", country)))

# aesthetics
clock <- readxl::read_excel("aesthetics/clock_scale.xlsx")
prettdir <- readxl::read_excel("aesthetics/clock_colors.xlsx")

#......................
# confirm that D9 and O6 are same except missing
#......................
check_cdmtdna <- seqinr::read.fasta("data/noclonefasta/unique_mtdna.fa")
check_cdmtdna <- check_cdmtdna[names(check_cdmtdna) %in% c("D9U3K", "O6Y4K")]
diffsites <- which(check_cdmtdna$D9U3K != check_cdmtdna$O6Y4K)
check_cdmtdna$O6Y4K[diffsites]
# as we can see from the viz, these are all just missing sites


#............................................................
# phylogeography
#...........................................................
# get hammings
hammdist <- ape::dist.dna(ape::as.DNAbin(mtdna_uniq), model = "N")
hammdist <- broom::tidy(hammdist) %>%
  dplyr::filter(item1 %in% c("D9U3K") |
                  item2 %in% c("D9U3K")) %>%
  dplyr::filter(!c(item1 %in% c("O6Y4K", "Q8J6O") |
                      item2 %in% c("O6Y4K", "Q8J6O")))
# expand grid
hammdist_expd <- tibble::tibble(item1 = hammdist$item2,
                                item2 = hammdist$item1,
                                distance = hammdist$distance)
# bring together
hammdist_sng <- dplyr::bind_rows(hammdist, hammdist_expd) %>%
  dplyr::filter(item1 == "D9U3K") %>%
  dplyr::rename(smpls = item2)

# make to clockwork
plotdat <- dplyr::left_join(hammdist_sng, mtdt, by = "smpls") %>%
  dplyr::select(c("item1", "smpls", "distance", "iid", "country")) %>%
  dplyr::left_join(., clock, by = "country") %>%
  dplyr::mutate(yclockposition = yscalar * yclockposition,
                distance = ifelse(distance == 0, 0.25, distance))

# add jitter
plotdat$ypos = plotdat$yclockposition + rnorm(nrow(plotdat), mean = 0, sd = 1)
# DRC data
drcdat <- tibble::tibble(distance = 0, ypos = 0)

# plot
plotObj <- plotdat %>%
  ggplot() +
  geom_point(data = drcdat, aes(x = distance, y = ypos),
             color = "#FF0018", shape = 8, size = 3) +
  geom_point(aes(x = distance, y = ypos, color = country, shape = country),
             alpha = 0.75, size = 2) +
  coord_polar(theta = "y") +
  xlab("Hamming's Distance") +
  scale_color_manual(name = "Country",
                     labels =c("BR","CN","CO","Ebro1944",
                               "ET","ID","IN","KH",
                               "LA", "LK",
                               "MG","MM",
                               "MX","MY","NHA",
                               "PE","PG",
                               "TH","VN"),
                     values = c(prettdir$hexcolor)) +
  scale_shape_manual(name = "Country",
                     labels =c("BR","CN","CO","Ebro1944",
                               "ET","ID","IN","KH",
                               "LA", "LK",
                               "MG","MM",
                               "MX","MY","NHA",
                               "PE","PG",
                               "TH","VN"),
                     values =  c(prettdir$shape)) +
  guides(color = guide_legend(nrow = 4, byrow = TRUE)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),
        axis.text.y = element_text(family = "Helvetica", size = 11, face = "italic"),
        legend.position = "bottom",
        legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.5, size = 12),
        legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10))


#......................
# bring together for manuscript
#......................
whitespace <- ggplot() +
  theme_void()

svglite::svglite("figures/clock_msn_figure.svg",
                 width = 10, height = 7)
cowplot::plot_grid(plotObj, whitespace,
                   nrow = 1, labels = c("(A)", "(B)"))
graphics.off()
