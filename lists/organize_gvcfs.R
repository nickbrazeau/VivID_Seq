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


