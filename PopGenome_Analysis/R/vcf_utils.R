#' @param vcflist list of s4 class from vcfR;
#' @details cbinds together without regard to pos
cbindvcflist <- function(vcflist){

  meta <- vcflist[[1]]@meta
  fix <- vcflist[[1]]@fix
  format <- vcflist[[1]]@gt[,1]

  vcflistgt <- purrr::map(vcflist, "gt")
  names(vcflistgt) <- NULL
  vcflistnofmt <- lapply(vcflistgt, function(x){return(x[,2:ncol(x)])})

  gt <- do.call("cbind.data.frame", vcflistnofmt)
  gt <- as.matrix( cbind(format, gt) )
  colnames(gt)[1] <- "FORMAT"

  # write out new vcfRobj
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)
  return(newvcfR)

}



#' @param vcfRobj s4 class from vcfR;
#' @details uses WSAF (via AD and DP) to force monoclonal calls

vcfR_force_hets_to_homo <- function(vcfRobj){

  # extract coverage and counts matrices
  coverage <- vcfR::extract.gt(vcfRobj, element = "DP", as.numeric = T)
  counts_raw <- vcfR::extract.gt(vcfRobj, element = "AD")
  ref <- masplit(counts_raw, record = 1, sort = FALSE, decreasing = FALSE)
  alt1 <- masplit(counts_raw, record = 2, sort = FALSE, decreasing = FALSE)

  # find het calls we need to lfitover
  hets <- vcfR::is.het(vcfR::extract.gt(vcfRobj))

  # get GT mats
  gtmat <- vcfR::extract.gt(vcfRobj)
  new.gtmat <- matrix(NA, nrow=nrow(gtmat), ncol=ncol(gtmat))

  opt <- c("ref", "alt1")

  for(j in 1:ncol(hets)){ # for every sample
    for(i in 1:nrow(hets)){ # for every loci
      if(hets[i,j]){
        refaf <- ref[i,j]/coverage[i,j]
        alt1af <- alt1[i,j]/coverage[i,j]

        swtch <- opt[ which( c(refaf, alt1af) ==
                               max(c(refaf, alt1af), na.rm = T)) ]

        if(length(swtch) > 1){ # multiple maxes, sample one
          swtch <- sample(swtch, size = 1)
        }


        switch(swtch,
               ref={
                 new.gtmat[i,j] <- "0/0"
               },
               alt1={
                 new.gtmat[i,j] <- "1/1"
               }
        )


      } else{
        new.gtmat[i,j] <- gtmat[i,j] # stays the same
      }
    }
  } #end nested for loop

  # liftover for vcf details
  gt <- ifelse(is.na(new.gtmat), NA,
               paste0(new.gtmat, ":", counts_raw, ":", coverage)) # loci as rows, smpls as columns

  # append format column and sample names
  gt <- cbind(FORMAT = "GT:AD:DP", gt)
  colnames(gt) <- colnames(vcfRobj@gt)

  # getFix
  fix <- vcfRobj@fix
  # get meta
  meta <- append(vcfRobj@meta, "##user=Forced this to homozygosity")

  # write out new vcfRobj
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)
  return(newvcfR)

}



#' @param vcfRobj s4 class from vcfR;
#' @details uses unique function to find GTs that are not identical across loci for all samples

vcfR_filter_unique_sites <- function(vcfRobj){
  gtmat <- vcfR::extract.gt(vcfRobj)
  segsites <- apply(gtmat, 1, function(x){return(length(unique(x)) != 1)})

  gt <- vcfRobj@gt[as.vector(segsites), ]
  fix <- vcfRobj@fix[as.vector(segsites), ]

  new.vcfRobj <- new("vcfR", meta = vcfRobj@meta, gt = gt, fix = fix)
  return(new.vcfRobj)
}
