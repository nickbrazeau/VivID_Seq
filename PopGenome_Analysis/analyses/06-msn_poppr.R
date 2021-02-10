####################################################################################
## Purpose: Minimum spanning network
##
## Notes:
####################################################################################
library(poppr)
library(igraph)
library(seqinr)
library(tidyverse)

#............................................................
# no clones fasta
#...........................................................
#uniq_mtdna <- ape::read.FASTA("data/noclonefasta/unique_mtdna.fa")
#uniq_mtdna <- uniq_mtdna[names(uniq_mtdna) != "O6Y4K"]
uniq_mtdna <- ape::read.FASTA("data/rep_smpls_noclones/rep_mtdna.fa")

#......................
# confirm that D9 and O6 are same except missing
#......................
check_cdmtdna <- seqinr::read.fasta("data/noclonefasta/unique_mtdna.fa")
check_cdmtdna <- check_cdmtdna[names(check_cdmtdna) %in% c("D9U3K", "O6Y4K")]
diffsites <- which(check_cdmtdna$D9U3K != check_cdmtdna$O6Y4K)
check_cdmtdna$O6Y4K[diffsites]
# as we can see from the viz, these are all just missing sites

#............................................................
# read in
#...........................................................
# aesthetics
prettdir <- readxl::read_excel("aesthetics/msn_colors.xlsx")

# metadat
mtdt <- readRDS("data/derived_data/mtdt.RDS")
mtdt <- mtdt %>%
  dplyr::mutate(country = ifelse(vividregion == "NHA", "NHA",
                                           country))
# ebro
mtdt <- mtdt %>%
  dplyr::mutate(country = ifelse(smpls == "Ebro1944", "Ebro1944", country))


# drop to samples were are considering
mtdt <- mtdt %>%
  dplyr::filter(smpls %in% names(uniq_mtdna))

#............................................................
# read in gendata
#...........................................................
genclone <- poppr::as.genclone(adegenet::DNAbin2genind(uniq_mtdna))
genclone
#......................
# add in metadata
#......................
adg_ord <- tibble::tibble(smpls = indNames(genclone))
adg_ord <- dplyr::left_join(adg_ord, mtdt, by = "smpls")
strata(genclone) <- adg_ord[,c("country", "vividregion")]
# set popopulations
setPop(genclone) <- ~country
genclone

#............................................................
# min span
#...........................................................
# distance w/ miss data
hammdist <- ape::dist.dna(uniq_mtdna,
                          model = "N")

# quick class sanity check
testdist <- provesti.dist(genclone)
identical(rownames(testdist), rownames(hammdist))
identical(colnames(testdist), colnames(hammdist))

#......................
# net
#......................
min_span_net <- poppr.msn(genclone, hammdist,
                          showplot = F,
                          include.ties = T)
# make new colors
color_order <- tibble::tibble(country = names(min_span_net$colors))
color_order <- color_order %>%
  dplyr::left_join(., prettdir, by = "country")

min_span_net$colors <- color_order$hexcolor


# save out
png("figures/minimum_spanning_network_rep_mtdna.png",
    width = 6, height = 6, units = "in", res  = 500)
plot_poppr_msn(genclone, min_span_net,
               inds = "none",
               mlg = F,
               wscale = T,
               nodelab = 999,
               nodescale = 100,
               quantiles = F,
               palette = min_span_net$colors,
               pop.leg = T,
               layout = layout_with_kk)
graphics.off()

