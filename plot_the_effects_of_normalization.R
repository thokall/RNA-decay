plot_the_effects_of_normalization <- function(DGEList = obj, specific.genes = c("sel", "aurA", "roX1"), random.genes = NA){
  # Plots the decay curve of some genes before and after normalization
  obj <- DGEList
  if(is.na(random.genes)){
    genes <- numeric()
    for(i in 1:length(specific.genes)){
      row <- which(obj$genes$symbol == specific.genes[i])
      if(length(row) == 0){
        stop(paste0("Gene name: " ,specific.genes[i], " not recognized."))
      }
      else{genes <- c(genes, which(obj$genes$symbol == specific.genes[i]))
      }
    }
  }
  if(!is.na(random.genes)){genes <- sample(1:nrow(obj),random.genes)
  }
  
  par(mfrow = c(length(genes),2))
  data <- obj$cpm_trimmed
  data  <- apply(data,2,as.numeric)
  norm_data <- obj$cpm_trimmed_normalized
  norm_data  <- apply(norm_data,2,as.numeric)
 
  ##create a time point vector that fits the main data
  t <- sapply(strsplit(colnames(data),"t"), function(x){x[!x ==""]})
  t <- as.numeric(t)
  
  for(i in 1:length(genes)){
    
    #adjust the time point vector to the specific gene
    realt <- suppressWarnings(t+data[i,]-data[i,])
    plot(realt,data[genes[i],], type = "b", main = paste("raw data:",obj$genes$symbol[genes[i]]), log = "y", xlab = "time")
    plot(realt,norm_data[genes[i],], type= "b", main = paste("normalized:",obj$genes$symbol[genes[i]]), log = "y", xlab = "time")
    
    
  
      

  }
  
  
   
  
}#closes function
  