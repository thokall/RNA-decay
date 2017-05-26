apply_expression_cutoff <- function(DGEList = obj, earlyCPM_cutoff = 2, CountCutoff = 10){
  ## Takes a DGEList object as an input and outputs the filtered DGEList object
  ## The samples must be named after the time point "t" followed by a number, eg "t0", "t1", "t2"
  ## filter on cpm of t0, t1, t2
  ## filter on sum of raw count
  ## recalculates library size
  ## prints to screen the number of genes filtered out and how many genes are left
  
  obj <- DGEList
  all_rows <- nrow(obj)
  keep <- rowSums(cpm(obj)) > CountCutoff
  obj <- obj[keep, , keep.lib.sizes = FALSE]
  countfilter <- nrow(obj)
  
  if(all(c("t0", "t1","t2") %in% colnames(cpm(obj)))){
    keep <- cpm(obj[,"t0"]) > earlyCPM_cutoff & cpm(obj[,"t1"]) > earlyCPM_cutoff & cpm(obj[,"t2"]) > earlyCPM_cutoff
    obj <- obj[keep, , keep.lib.sizes = FALSE]
    cpmfilter <- nrow(obj)
  }
  else if(all(c("t0", "t1") %in% colnames(cpm(obj)))){
    keep <- cpm(obj[,"t0"]) > earlyCPM_cutoff & cpm(obj[,"t1"]) > earlyCPM_cutoff 
    obj <- obj[keep, , keep.lib.sizes = FALSE]
    cpmfilter <- nrow(obj)
  }
  else if(all(c("t0","t2") %in% colnames(cpm(obj)))){
    keep <- cpm(obj[,"t0"]) > earlyCPM_cutoff & cpm(obj[,"t2"]) > earlyCPM_cutoff
    obj <- obj[keep, , keep.lib.sizes = FALSE]
    cpmfilter <- nrow(obj)
  }
  else if(all(c("t1","t2") %in% colnames(cpm(obj)))){
    keep <- cpm(obj[,"t1"]) > earlyCPM_cutoff & cpm(obj[,"t2"]) > earlyCPM_cutoff
    obj <- obj[keep, , keep.lib.sizes = FALSE]
    cpmfilter <- nrow(obj)
  }
  else{stop("There are not enough early time points ... exiting ...")}
  
  cat("There were ", all_rows, " genes in the input. \n", "After filtering for sum of row counts (", CountCutoff, "), there are ",countfilter, " left. \n", "After filtering for CPM at t=0 (", earlyCPM_cutoff,"), there are ", cpmfilter, " genes left. \n", sep="")
  return(obj)
}