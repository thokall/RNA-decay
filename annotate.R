annotate <- function(DGEList = obj, organism = c("dmel","other")){
  ## requires access to internet
  ## Takes a DGEList object as an input and outputs the annotated DGEList object
  ## Annotates with chromosome name, symbol, start position, stop position
  ## only available for drosophila melanogaster at the moment
  if(organism == "dmel"){
    library("biomaRt")
    #mart <- useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
    mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl", host="www.ensembl.org")
    annots <- getBM(attributes=c("flybase_gene_id", "chromosome_name","external_gene_name","start_position","end_position"), filters="flybase_gene_id", values= obj$genes, mart=mart)
    m <- match(obj$genes$genes,annots$flybase_gene_id)
    obj$genes$chromosome_name <- annots$chromosome_name[m]
    obj$genes$symbol <- annots$external_gene_name[m]
    obj$genes$start_position <- annots$start_position[m]
    obj$genes$end_position <- annots$end_position[m]
  }
  else{print("This organism is not supported")}
  obj
}


