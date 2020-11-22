vividbams <- list.files(path = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/wgs_pe_improved_ViVIDSmpls/aln/merged",
                        full.names = T, pattern = ".bam")
vividbams <- vividbams[!grepl(".bai", vividbams)]

# symlinked ebro under pe for future ease
# sebams <- list.files(path = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/wgs_se_improved_global/aln/merged",
#                         full.names = T, pattern = ".bam")
# sebams <- sebams[!grepl(".bai", sebams)]

pebams <- list.files(path = "~/Documents/MountPoints/mountedScratchLL/Projects/VivID_Seq/wgs_pe_improved_global/aln/merged",
                     full.names = T, pattern = ".bam")
pebams <- pebams[!grepl(".bai", pebams)]

all_bams <- tibble::enframe( c(vividbams, sebams, pebams), name = NULL )

# fix local to remote
all_bams$value <- gsub("/Users/nickbrazeau/Documents/MountPoints/mountedScratchLL/",
                       "/pine/scr/n/f/nfb/", all_bams$value)



readr::write_tsv(x = all_bams,
                 path = "../lists/all_bams.list",
                 col_names = F)


