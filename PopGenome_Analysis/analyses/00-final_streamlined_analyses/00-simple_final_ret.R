## .................................................................................
## Purpose: Simple Script to bring together final results
##          to include for phylogenetic analysis
##
## Notes: This is the only file that matters in this analysis repo.
##        The other RMD scripts are depreciated but kept along for
##        those interested in the additional data
## .................................................................................
library(tidyverse)
library(PopGenome)
library(Biostrings)

#............................................................
###### Section 1: Quick stats #####
#...........................................................
#......................
# in data
#......................
mtdt <- readRDS("data/derived_data/mtdt.RDS")
# long vcf
vcf.snp.ss.long_gt <- readRDS("data/derived_data/final_vcf_snps_long.RDS")
vcf.snp.ss.long_gt <- vcf.snp.ss.long_gt %>%
  dplyr::left_join(., y = mtdt, by = c("iid", "country", "vividregion")) %>%
  dplyr::filter(passed == "Y") %>%
  dplyr::mutate(POS = factor(POS, levels = c(1:5989), ordered = T))

# all dna
mtdna_all <- Biostrings::readDNAStringSet(filepath = "data/fasta/mtdna_anc.fa")
# no clones
mtdna_uniq <- Biostrings::readDNAStringSet(filepath = "data/noclonefasta/unique_mtdna.fa")


#......................
# basic popgen statistics
#......................
mtdna_uniq
mtdna_all
# get seg sites
seqmat <-  tibble::as_tibble( t(as.matrix(mtdna_uniq) ) )
seqmat <- cbind.data.frame(POS = 1:nrow(seqmat), seqmat)
segsites <- apply(seqmat[,2:ncol(seqmat)], 1, function(x){
  ret <- unique(x)
  ret <- ret[ret != "N"]
  return(ret)
})
sum(sapply(segsites, length) != 1)

# define countries
mtdt <- mtdt %>%
  dplyr::mutate(effcountry = ifelse(vividregion == "Lab", "Lab",
                                    ifelse(vividregion == "NHA", "NHA",
                                           country)))
mtdt_noclones <- mtdt %>%
  dplyr::filter(smpls %in% names(mtdna_uniq))
table(mtdt_noclones$effcountry)

# read in and set up popgenome object
noclones_popgen <- PopGenome::readData(path = "data/noclonefasta/",
                                       format = "fasta",
                                       include.unknown = T,
                                       progress_bar_switch=T,
                                       FAST=F,
                                       big.data=F)
rgn <- split(mtdt_noclones$smpls, factor(mtdt_noclones$effcountry))
noclones_popgen <- PopGenome::set.populations(noclones_popgen,
                                              rgn)
# set Cynmologi as outgroup
noclones_popgen <- PopGenome::set.outgroup(noclones_popgen,
                                           c("ERS001838")) #PcyM-strain
get.sum.data(noclones_popgen)

#......................
# w/in diversity
#......................
# must be calculated first before diversity and neutrality
noclones_popgen <- PopGenome::F_ST.stats(noclones_popgen)
noclones_popgen <- PopGenome::diversity.stats.between(noclones_popgen)
noclones_popgen <- PopGenome::neutrality.stats(noclones_popgen)
ret <- PopGenome::get.diversity(noclones_popgen, between = F)
popn <- gsub(" ", "", rownames(ret))
tbout <- lapply(ret, as.data.frame) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(popn = names(noclones_popgen@populations)) %>%
  dplyr::select(c("popn", dplyr::everything()))

# write out
dir.create("tables")
readr::write_csv(tbout, file = "table/basic_popgen_summstat.csv")

#............................................................
# basic hammings distance
#...........................................................
hammdist <- ape::dist.dna(ape::as.DNAbin(mtdna_uniq), model = "N")
zerodist <- hammdist %>%
  broom::tidy(.) %>%
  dplyr::filter(distance == 0) %>%
  dplyr::filter(item1 %in% c("D9U3K", "O6Y4K", "Q8J6O") |
                  item2 %in% c("D9U3K", "O6Y4K", "Q8J6O"))

zerodist <- unique(c(as.character(zerodist$item1), as.character(zerodist$item2)))
zerodist <- zerodist[!zerodist %in% c("D9U3K", "O6Y4K", "Q8J6O")]
# get smpls
zerodist_smpls <- mtdt %>%
  dplyr::filter(smpls %in% zerodist)

#......................
# nha and ebro dist
#......................
nhasmpls <- mtdt %>%
  dplyr::filter(vividregion == "NHA") %>%
  dplyr::pull("smpls")
hammdist %>%
  broom::tidy(.) %>%
  dplyr::filter(item1 %in% c("D9U3K", "O6Y4K", "Q8J6O") |
                  item2 %in% c("D9U3K", "O6Y4K", "Q8J6O")) %>%
  dplyr::filter(item1 %in% c(nhasmpls, "Ebro1944") |
                  item2 %in% c(nhasmpls, "Ebro1944"))

#......................
# all
#......................
alldrcpairs <- hammdist %>%
  broom::tidy(.) %>%
  dplyr::filter(item1 %in% c("D9U3K", "O6Y4K", "Q8J6O") |
                  item2 %in% c("D9U3K", "O6Y4K", "Q8J6O")) %>%
  dplyr::filter(!(item1 == "ERS001838" |
                    item2 == "ERS001838")) # drop PyCm
summary(alldrcpairs)

#............................................................
###### Section 2: GT Viz #####
#...........................................................
# find seg sites from fasta
mtdna_uniq_nocym <- mtdna_uniq[names(mtdna_uniq) != "ERS001838"]
seqmat <-  tibble::as_tibble( t(as.matrix(mtdna_uniq_nocym) ) )
seqmat <- cbind.data.frame(POS = 1:nrow(seqmat), seqmat)
segsites <- apply(seqmat[,2:ncol(seqmat)], 1, function(x){ return( length(unique(x)) > 1 ) })

# take to long
seqmat.long <- seqmat %>%
  dplyr::filter(segsites) %>%
  tidyr::gather(., key = "smpls", value = "GT", 2:ncol(.)) %>%
  dplyr::mutate(GT = ifelse(GT == "N", NA, GT)) %>%
  dplyr::mutate(GT = factor(GT)) %>%
  dplyr::left_join(., mtdt, by = "smpls") %>%
  dplyr::filter(hapuid == "Y") %>%
  dplyr::mutate(POS = factor(POS, levels = c(1:5989), ordered = T))


# make plot
gtplot <- seqmat.long %>%
  ggplot() +
  geom_tile(aes(x = smpls, y = POS, fill = GT)) +
  scale_fill_manual("Alleles", values = c("#4daf4a", "#377eb8", "#ffff33", "#e41a1c")) +
  ylab("POS") +
  facet_wrap(~effcountry, scales = "free_x", nrow = 1, strip.position = "top") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_blank())

#............................................................
###### Section 3: Phylogen #####
#...........................................................
set.seed(48)
library(ape)
library(phangorn)
mtdna.ape <- ape::as.DNAbin(mtdna_uniq)
mtdna.dist.JC <- ape::dist.dna(mtdna.ape, model="JC69")
# Phangorn
mtdna.phangorn <- phangorn::as.phyDat(mtdna.ape)
tree.init <- ape::nj(mtdna.dist.JC)

# JC Model
JCfit.init <- phangorn::pml(tree.init, mtdna.phangorn)
fitJC <- phangorn::optim.pml(JCfit.init, model = "JC",
                             rearrangement = "stochastic",
                             optNni=TRUE, optBf=TRUE, optQ=TRUE, optGamma=TRUE)

# GTR Model
fitGTR.init <- phangorn::pml(tree.init, mtdna.phangorn, k=4)
fitGTR <- phangorn::optim.pml(fitGTR.init, model = "GTR",
                              rearrangement = "stochastic",
                              optNni=TRUE, optBf=TRUE, optQ=TRUE, optGamma=TRUE)

# compare two models
out <- cbind.data.frame(model = c("JC69", "GTR"),
                        as.data.frame(anova(fitJC, fitGTR)),
                        AIC = c(AIC(fitJC), AIC(fitGTR)))
knitr::kable(out)

# Bootstrap GTR which is better model
iters <- 5000
GTRbs <- phangorn::bootstrap.pml(fitGTR, bs=iters,
                                 optNni=TRUE,
                                 rearrangement = "stochastic")

#...........................................................
# plot tree
#...........................................................
library(ggtree)
#......................
# organize bootstrap
#......................
bstree <- plotBS(midpoint(fitGTR$tree), GTRbs, p = 0, type="p")

node_labels_in_edge <- bstree$node.label[bstree$edge[,1]-Ntip(bstree)]
tips_nodes <- bstree$edge[,2]

#......................
# make prettty tree
#......................
prettdir <- readxl::read_excel("aesthetics/tree_color_decode.xlsx")
tidytree <-  treeio::as.treedata(bstree)

# add in boots
bp2 <- tibble(node=1:Nnode(fitGTR$tree) + Ntip(fitGTR$tree),
              bootstrap = prop.clades(fitGTR$tree, GTRbs, rooted = F))
tidytree <- full_join(tidytree, bp2, by="node")

# tidy and make tree
treemtdt <- mtdt %>%
  dplyr::mutate(effcountry = ifelse(smpls == "ERS001838", "pcynomolgi", effcountry)) %>%
  dplyr::rename(label = smpls) %>%
  dplyr::filter(hapuid == "Y") %>%
  dplyr::select(c("label", "effcountry")) %>%
  dplyr::mutate(treegroup = ifelse(label == "Ebro1944", "Ebro1944",
                                   ifelse(label == "ERS001838", "Pcynomolgi",
                                          ifelse(effcountry == "CD", "DRC", effcountry))))

tidytree <- full_join(tidytree, treemtdt, by = 'label')
tidytree@extraInfo$treegroup[is.na(tidytree@extraInfo$treegroup)] <- "zinternal"
tidytree@extraInfo$treegroup[is.na(tidytree@extraInfo$treegroup)] <- "zinternal"

# extra manip to trick ggplot into combining legends
tidytree@extraInfo <- bind_rows(tidytree@extraInfo,
                                tibble::tibble(node = 999,
                                               effcountry = "zinternal",
                                               treegroup = "zinternal"))

# make plot
bstree_plotObj <- ggtree(tidytree,
                         branch.length="none",
                         layout="circular",
                         aes(color = treegroup)) +
  geom_tippoint(aes(color =  treegroup, shape = treegroup),
                show.legend = T,
                size = 1.7) +
  geom_nodelab() +
  geom_point2(aes(subset=!isTip,
                  color=cut(bootstrap, c(0, 400, 500, 600, 700, Inf))),
              shape=21, size=3, color = "transparent") +
  scale_fill_manual(values = rev(c("#f0f0f0", "#d9d9d9", "#969696", "#737373", "#252525")),
                    guide='legend',
                    name='Bootstrap Percentage (BP)',
                    breaks=c('(700,1e+03]', '(600,700]', '(500,600]', '(400,500]', '(0,400]'),
                    labels=expression(BP>=70,
                                      60 <= BP * " < 70",
                                      50 <= BP * " < 60",
                                      40 <= BP * " < 50",
                                      BP < 40)) +
  scale_color_manual(name = "Country",
                     labels =c("Brazil","China","Colombia","DRC","Ebro1944",
                               "Ethiopia","Indonesia","India","Cambodia",
                               "Laos","Lab","Sri Lanka",
                               "Madagascar","Myanmar",
                               "Mexico","Malaysia","NHA",
                               "Pcynomolgi","Peru","Papua New Guinea",
                               "Thailand","Vietnam", "zinternal"),
                     values = c(prettdir$hexcolor, "#000000")) +
  scale_shape_manual(name = "Country",
                     labels =c("Brazil","China","Colombia","DRC","Ebro1944",
                               "Ethiopia","Indonesia","India","Cambodia",
                               "Laos","Lab","Sri Lanka",
                               "Madagascar","Myanmar",
                               "Mexico","Malaysia","NHA",
                               "Pcynomolgi","Peru","Papua New Guinea",
                               "Thailand","Vietnam", "zinternal"),
                     values =  c(prettdir$shape, 0)) +
  theme(legend.title = element_text(size = 12,
                                    hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
  )


svglite::svglite("figures/no_clones_tree.svg",
                 width = 11, height = 8)
bstree_plotObj
graphics.off()

#......................
# bs values for reference
#......................
tree_bsvals <- ggtree(tidytree,
                      branch.length="none",
                      aes(color = treegroup)) +
  geom_tippoint(aes(color =  treegroup, shape = treegroup), show.legend = F, size = 0.8) +
  geom_nodelab(aes(subset=!isTip,
                   x=branch, label=bs)) +
  scale_color_manual(name = "Country",
                     labels =c("Brazil","China","Colombia","DRC","Ebro1944","Ethiopia","Indonesia" ,"India","Cambodia","Laos","Lab","Sri Lanka",
                               "Madagascar","Myanmar","Mexico","Malaysia","NHA","Pcynomolgi","Peru","Papua New Guinea","Thailand","Vietnam", "zinternal"),
                     values = c(prettdir$hexcolor, "#000000")) +
  scale_shape_manual(name = "Country",
                     labels = c("Brazil","China","Colombia","DRC","Ebro1944","Ethiopia","Indonesia" ,"India","Cambodia","Laos","Lab ","Sri Lanka",
                                "Madagascar","Myanmar","Mexico","Malaysia","NHA","Pcynomolgi","Peru","Papua New Guinea","Thailand","Vietnam", "zinternal"),
                     values =  c(prettdir$shape, 0)) +
  theme_tree() +
  theme(legend.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold"))

tree_bsvals


jpeg("figures/no_clones_bootstrap_values_tree.jpg", width = 11, height = 8, units = "in", res = 600)
tree_bsvals
graphics.off()
