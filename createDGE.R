createDGE <- function(countMatrix = NA, annotationFile = "annotation.gtf", sampleInfo = "targets.txt", strand = 0, paired = c("YES", "NO"), nthreads = 1) {
  ## Create a DGE list by counting reads on bam files that are
  ## listed in the sampleInfo dataframe. With no arguments specified
  ## the function reads the file targets.txt and  annotation.gtf in the
  ## current directory and assumes that the data is unstranded.
  ##
  ## Args:
  ## countMatrix: Character vector length 1. Name and path to available count matrix.
  ## The file should have column header and the first column should be gene names.
  ## annotationFile: Character vector length 1. The complete path to the gtf file used
  ## for analysis. Annotation.gtf in the current folder will be used as default
  ## sampleInfo: Character vector length 1. Name of the file containing
  ## information on the bam files to be analysed.
  ## This file needs to have a column named "sample_name" containing
  ## filenames and path to all input files
  ## strand: Numeric vector of length 1. Is the data stranded or not.
  ## 0 = unstranded,
  ## 1 = stranded with first read in direction of annotation,
  ## 2 = stranded with first read opposite of annotation (Typical for Illumina)
  ## paired: Character vector length 1. Is the data paired end (YES), or not (NO)
  ## nthreads: Numeric vector length 1. Number of threads used for counting reads.
  wd <- getwd()
  if (is.na(countMatrix)) {
    pe = match.arg(paired, c("YES", "NO"))
    sampleInfo <- read.table(sampleInfo, header = TRUE)
    bamFiles = sampleInfo$sample_name
    fc <- featureCounts(files = bamFiles, annot.ext = annotationFile,
                        isPaired = pe, nthreads = nthreads,
                        isGTFAnnotationFile = TRUE, strandSpecific = strand)
    dge <- DGEList(counts = fc$counts, genes = fc.annotation)
    dge
  } else {
    cm <- read.table(paste0(wd,"/",countMatrix), header = TRUE)
    dge <- DGEList(counts = cm[,-1], genes = cm[,1])
    dge
  }
}