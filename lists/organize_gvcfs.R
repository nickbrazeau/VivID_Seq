###############################################
####        Read in mtdt & gvcfs           ####
###############################################
smpls <- readxl::read_excel("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/vivid_seq_public_NGS.xlsx")

all_gvcfs <- tibble::enframe( list.files(path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/vcfs_gatk_joint_raw",
                                         pattern = ".g.vcf",
                                         full.names = T),
                              name = NULL )

# fix local to remote
all_gvcfs$value <- gsub("/Users/nickbrazeau/Documents/MountPoints/mountedScratchLL/",
                       "/pine/scr/n/f/nfb/", all_gvcfs$value)


###############################################
#### Read in Metadata & Write Lab strains  ####
###############################################
all_gvcfs$basenames <- gsub(".bam", "", basename(all_gvcfs$value))

lab_strains <- all_gvcfs$value[ all_gvcfs$basenames %in% smpls$acc[smpls$host == "Lab_strain"]]
lab_strains <- tibble::enframe( lab_strains, name = NULL )

readr::write_tsv(x = lab_strains,
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/lists/lab_strains.gvcfs.list",
                 col_names = F)



###############################################
#### Get Hybrid Capture Sample Comparison  ####
###############################################
# KP053-HYB,   KP063-HYB,   OM032-HYB,  OM092-HYB
# SRS1061008,  SRS1061007,  SRS1061044, SRS1061081
hybridcompar <- c("KP053-HYB", "KP063-HYB", "OM032-HYB", "OM092-HYB",
                  "SRS1061008", "SRS1061007", "SRS1061044", "SRS1061081")


hybridcompar <- all_gvcfs$value[ all_gvcfs$basenames %in% hybridcompar ]
hybridcompar <- tibble::enframe( hybridcompar, name = NULL )

readr::write_tsv(x = hybridcompar,
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/lists/hybridcompar.gvcfs.list",
                 col_names = F)

