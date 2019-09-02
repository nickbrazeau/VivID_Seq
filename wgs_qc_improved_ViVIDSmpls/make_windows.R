fai <- readr::read_tsv("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Seq/wgs_qc_improved/chroms.bed", col_names = F)


make_chrom_bed_windows <- function(chrombed, windowsize = 5e3){
  chrombed.wind <- split(chrombed, 1:nrow(chrombed)) %>% 
    lapply(., function(x){
      
      start <- 0
      pos <- windowsize
      ret <- NULL
      
      while(pos < x[,3]){
        ret <- rbind.data.frame(ret, 
                                data.frame(chrom = unlist( unique(x[,1]) ), 
                                           start = start,
                                           end = pos
                                ))
        start <- pos
        pos <- pos + windowsize
        
      }
      ret <- rbind.data.frame(ret, 
                              data.frame(chrom = unlist( unique(x[,1]) ), 
                                         start = start,
                                         end = unlist(x[,3]) )
      ) # end of chrom
      
      ret$chrom <- paste0(ret$chrom, "_chunk", 1:nrow(ret))
      
    return(ret)
        
    }) %>% 
    dplyr::bind_rows()
  
  return(chrombed.wind)
}


chrombed.wind <- make_chrom_bed_windows(fai, windowsize = 1e3)
