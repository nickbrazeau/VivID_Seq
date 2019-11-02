library(tidyverse)
###############################################
####             Read in mtdt              ####
###############################################
smpls <- readxl::read_excel("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/vivid_seq_public_NGS.xlsx")

###############################################
####              Read in QC               ####
###############################################
# criterion that 95% of mitogenome must be "callable"
callable <- round(5989*0.95, 0)
read_callable_loci <- function(path){
  ret <- readr::read_tsv(path, col_names = T)
  ret$smpl <- basename(path)
  ret$smpl <- gsub(".callable_summary.txt", "", ret$smpl)

  ret <- ret %>%
    dplyr::filter(state == "CALLABLE")

  return(ret)

}

globalqc <- list.files(path = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/wgs_qc_improved_global/qc/",
                       pattern = ".callable_summary.txt", full.names = T)
globalqc <- globalqc[!grepl("all.callable_summary.txt", globalqc)]
vividqc <- list.files(path = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/wgs_qc_improved_ViVIDSmpls/qc/",
                       pattern = ".callable_summary.txt", full.names = T)
qcpaths <- c(globalqc, vividqc)

qc.results <- lapply(qcpaths, read_callable_loci) %>%
  dplyr::bind_rows()

qc.passed <- qc.results %>%
  dplyr::filter(nBases >= callable) %>%
  dplyr::select(c("smpl")) %>%
  unlist(.)




###############################################
####     Read in paths for gvcfs           ####
###############################################
all_gvcfs <- tibble::enframe( list.files(path = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/vcfs_gatk_joint_raw/chunks/all/",
                                         pattern = ".g.vcf.gz$",
                                         full.names = T),
                              name = NULL )

# fix local to remote
all_gvcfs$value <- gsub("/Users/nickbrazeau/Documents/MountPoints/mountedScratchLL/",
                       "/pine/scr/n/f/nfb/", all_gvcfs$value)

all_gvcfs$basenames <- gsub(".g.vcf.gz", "", basename(all_gvcfs$value))



###############################################
####         Write Lab strains             ####
###############################################
# note, I know all lab strains pass callable threshold
lab_strains <- all_gvcfs$value[ all_gvcfs$basenames %in% smpls$acc[smpls$host == "Lab_strain"]]
lab_strains <- tibble::enframe( lab_strains, name = NULL )

readr::write_tsv(x = lab_strains,
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/lists/lab_strains.gvcfs.list",
                 col_names = F)


###############################################
####      Write Filtered Samples           ####
###############################################
passed_smpls <- all_gvcfs$value[ !all_gvcfs$basenames %in% smpls$acc[smpls$host == "Lab_strain"] &
                                   all_gvcfs$basenames %in% qc.passed
                                   ]
passed_smpls <- tibble::enframe( passed_smpls, name = NULL )

readr::write_tsv(x = passed_smpls,
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/lists/passed_smpls.gvcfs.list",
                 col_names = F)



