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
# Read in lab strains
#........................................................
labstrains <- readr::read_tsv("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/lists/lab_strains.txt",
                              col_names = F)
labstrains <- gsub(".bam", "", basename(labstrains$X1))

#........................................................
# Read in vars
#........................................................
vcf <- vcfR::read.vcfR("~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/vcfs_gatk_joint_raw/all_raw.vcf.gz")
labstrain.vcf <- vcfRmanip::select_samples(vcfRobject = vcf, smplvctr = labstrains)


labstrain.vcf.ss <- vcfRmanip::vcfR2segsites_gt(labstrain.vcf)



labstrain.vcf.dp <- vcfR::extract.gt(labstrain.vcf, element = "DP")
labstrain.vcf.dp <- apply(labstrain.vcf.dp, 2, function(x){return(as.numeric(x))})


labstrain.vcf.dp %>%
  as_tibble(.) %>%
  tidyr::gather(., key = "smpl", value = "DP") %>%
  ggplot() +
  geom_boxplot(aes(x=factor(smpl), y = DP))


