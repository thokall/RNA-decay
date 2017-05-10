# This would be an example not using functions. In the function I have
# not yet added gtf files as attribute, and it is more cosmetics so it
# does not really alter the use.

sampleInfo <- read.table("pilot_conf.txt", header = TRUE)
attr(sampleInfo, "GTF") <- "/proj/b2015023/private/genomes/dmel-all-r6.12.gtf" # add gtf file info as attribute
fc <- featureCounts(files = sampleInfo$sample_name, annot.ext = attr(sampleInfo, "GTF"), isPairedEnd = TRUE, nthreads = 2, isGTFAnnotationFile = TRUE, strandSpecific = TRUE)
save.image(file = "FC.rdata")

# The function needs the packages edgeR, limma, Rsubread. For proper use one can try and capture if these are installed or not, but for know one has to manually do so.

if(require(Rsubread)) {
    message('Rsubread loaded correctly')
} else {
    install.packages('Rsubread')
}

