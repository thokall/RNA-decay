apply_expression_cutoff <- function(DGE = obj, earlyCut = 2, expressionCut = 10, earlyCols = 3) {
  ## Takes a DGEList object as an input and outputs the filtered
  ## DGEList object, with an extra slot named summaryFilter that
  ## contains the same information as the summary shown on screen
  ## after running the script. It first filters based on the total cpm
  ## for a gene over all time points. It can also filter on the early
  ## time points independently as genes with very low expression at
  ## start of the experiment are not informative on decay patterns
  ## over time.
    obj <- DGE
    all_rows <- nrow(obj)
    keep <- rowSums(cpm(obj)) > expressionCut
    obj <- obj[keep, , keep.lib.sizes = FALSE]
    countfilter <- nrow(obj)

    if(earlyCols < 2) {
        stop("There are not enough early time points ... exiting ...")
    } else {
        keep <- as.vector(apply(cpm(obj)[,1:earlyCols] > earlyCut, 1, all))
        obj <- obj[keep, , keep.lib.sizes = FALSE]
        cpmfilt <- nrow(obj)
    }

    summaryFilter <- data.frame(row.names = c("Genes at start:", "Cutoff for row sums:", "Genes after row sum filter:", "Number of early time points:", "Cutoff for early time points:","Genes after early cutoff:"), c(all_rows, expressionCut, countfilter, earlyCols, earlyCut, cpmfilt))
    colnames(summaryFilter) <- NULL
    cat("There were ", all_rows, " genes in the input. \n", "After filtering for sum of row counts (", expressionCut, "), there are ",countfilter, " left. \n", "After filtering for CPM at t=0 (", earlyCut,"), there are ", cpmfilt, " genes left. \n", sep="")
    obj$summaryFilter <- summaryFilter
    obj
}


