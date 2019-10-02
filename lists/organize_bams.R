vividbams <- list.files(path = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/wgs_pe_improved_ViVIDSmpls/aln/merged",
                        full.names = T, pattern = ".bam")
vividbams <- vividbams[!grepl(".bai", vividbams)]

sebams <- list.files(path = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/wgs_se_improved_global/aln/merged",
                        full.names = T, pattern = ".bam")
sebams <- sebams[!grepl(".bai", sebams)]

pebams <- list.files(path = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/wgs_pe_improved_global/aln/merged",
                     full.names = T, pattern = ".bam")
pebams <- pebams[!grepl(".bai", pebams)]

all_bams <- tibble::enframe( c(vividbams, sebams, pebams), name = NULL )

# fix local to remote
all_bams$value <- gsub("/Users/nickbrazeau/Documents/MountPoints/mountedScratchLL/",
                       "/pine/scr/n/f/nfb/", all_bams$value)



readr::write_tsv(x = all_bams,
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/lists/all_bams.list",
                 col_names = F)


###############################################
#### Read in Metadata & Write Lab strains  ####
###############################################
smpls <- readxl::read_excel("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/scrape_pubseqs/vivid_seq_public_NGS.xlsx")
all_bams$basenames <- gsub(".bam", "", basename(all_bams$value))

lab_strains <- all_bams$value[ all_bams$basenames %in% smpls$acc[smpls$host == "Lab_strain"]]
lab_strains <- tibble::enframe( lab_strains, name = NULL )

readr::write_tsv(x = lab_strains,
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/lists/lab_strains.list",
                 col_names = F)



###############################################
#### Get Hybrid Capture Sample Comparison  ####
###############################################
# KP053-HYB,   KP063-HYB,   OM032-HYB,  OM092-HYB
# SRS1061008,  SRS1061007,  SRS1061044, SRS1061081
hybridcompar <- c("KP053-HYB", "KP063-HYB", "OM032-HYB", "OM092-HYB",
                  "SRS1061008", "SRS1061007", "SRS1061044", "SRS1061081")


hybridcompar <- all_bams$value[ all_bams$basenames %in% hybridcompar ]
hybridcompar <- tibble::enframe( hybridcompar, name = NULL )

readr::write_tsv(x = hybridcompar,
                 path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/lists/hybridcompar.list",
                 col_names = F)

