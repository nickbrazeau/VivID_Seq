#################################################################################
#### Purpose of this script is subset to the lab strain vcfs
#### and clean for "gold standard" SNPs.
#### Given the depth at which these lab strains were sequenced,
#### we should believe the GATK calls without many hard filters
#################################################################################
library(tidyverse)
library(vcfR)
library(vcfRmanip)


#........................................................
# Read in vars
#........................................................
labstrain.vcf <- vcfR::read.vcfR("~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/vcfs_gatk_lab_strains/Pvlabstrains_filtered.vcf.gz")

# most are biallelic
labstrain.vcf.bi <- labstrain.vcf[is.biallelic(labstrain.vcf), ]
labstrain.vcf.bi
# all bi are segregating sites
labstrain.vcf.ss <- vcfRmanip::vcfR2segsites_gt(labstrain.vcf.bi)
labstrain.vcf.ss

# reasonable depth for var
# but some labstrains lacking
labstrain.vcf.dp <- vcfR::extract.gt(labstrain.vcf, element = "DP")
labstrain.vcf.dp <- apply(labstrain.vcf.dp, 2, function(x){return(as.numeric(x))})
labstrain.vcf.dp %>%
  as_tibble(.) %>%
  tidyr::gather(., key = "smpl", value = "DP") %>%
  ggplot() +
  geom_boxplot(aes(x=factor(smpl), y = DP)) +
  theme(axis.text.x = element_text(angle = 90))

summary(labstrain.vcf.dp)



# fitler vcfR
chromR <- vcfR::create.chromR(labstrain.vcf, "PvP01_MIT_v1")
plot(chromR)

chromR <- proc.chromR(chromR, verbose=TRUE)
plot(chromR)


#........................................................
# Viz and Critique SNPs
#........................................................
labstrain.vcf.gt <- vcfR::extract_gt_tidy(labstrain.vcf)

labstrain.vcf.gt %>%
  dplyr::group_by(Indiv) %>%
  dplyr::mutate(POS = vcfR::getPOS(labstrain.vcf),
                POS = factor(POS, levels = c(1:5898), ordered = T)) %>%
  ggplot() +
  geom_tile(aes(x = Indiv, y = POS, fill = gt_GT)) +
  theme(axis.text.x = element_text(angle = 90))





