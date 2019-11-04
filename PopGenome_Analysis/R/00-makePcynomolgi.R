#----------------------------------------------------------------------------------------------------
# Purpose of this script is to make a P.cynomologi Conensus Script
#----------------------------------------------------------------------------------------------------
library(seqinr)
library(vcfR)
library(vcfRmanip)
library(tidyverse)
source("R/vcf_utils.R")
set.seed(48)

#........................................................
# Read in filtered Pcynomologi
#........................................................
pc <- vcfR::read.vcfR("~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/pcynomolgi_vars/passed.joint.vcf.gz")
pc <- vcfR::extract.indels(pc, return.indels = F)

#........................................................
# Quick View
#........................................................
pcgtmat <- cbind.data.frame(POS = vcfR::getPOS(pc), vcfR::extract.gt(pc, element = "GT"))
pcgtmat %>%
  tidyr::gather(., key = "smpl", value = "GT", 2:ncol(.)) %>%
  dplyr::mutate(loci = paste0("mt", POS)) %>%
  ggplot() +
  geom_tile(aes(x=factor( loci ), y=factor(GT))) +
  facet_wrap(~smpl) +
  theme(axis.text.x = element_text(angle = 90))
# very consistent


#........................................................
# No Hets when Just look at SNPs
# Can just talk Alt and mutate PvP01 backbone from there
#........................................................
POS <- vcfR::getPOS(pc)
ALTS <- vcfR::getALT(pc)


#........................................................
# Read in PvP01 backbone and mutate
#........................................................
pvp01 <- seqinr::read.fasta("~/Documents/MountPoints/mountIDEEL/resources/genomes/Pvivax/genomes/PvP01.fasta",
                            forceDNAtolower = F)
pcfasta <- pvp01 <- pvp01$PvP01_MIT_v1

pcfasta[POS] <- ALTS

#........................................................
# Write out
#........................................................
seqinr::write.fasta(sequences = pcfasta, names = "Pcynomolgi",
                    file.out = "data/fasta/Pcynomolgi.fasta")


