calculate_half_life <- function(DGEList = obj){
  
  data <- DGEList$cpm_trimmed_normalized
  ##create a time point vector that fits the main data
  t <- sapply(strsplit(colnames(data),"t"), function(x){x[!x ==""]})
  
  results <- as.data.frame(data[0,])
  
  
  
  for(i in 1:nrow(data)){
    # adjust the length of the time point vector to the actual data
    t <- as.numeric(t)
    realt <- suppressWarnings(t+data[i,]-data[i,])
    
    #logtranform data
    log2data <- apply(data,2,log2)
    
    #fit_nls
    values <- data[i,]
    fit_nls = tryCatch(nls(values ~ a*exp(-b*realt), data=list(realt,values), start=list(a= data[i,1],b=0.1)), error=function(e){paste("nls fitting failed row:" ,i)})
    if(grepl("fitting failed",fit_nls)){
      b <- NA
      hl_nls <- NA
    }
      
      else{
        b <- coef(fit_nls)["b"]
        hl_nls <- log(2)/b
      }
    
    #fit_lm
    values <- log2data[i,]
    fit_lm <- tryCatch(lm(formula =  values ~ realt), error=function(e){paste("lm fitting failed row:" ,i)})
    if(grepl("fitting failed",fit_lm)){
      hl_lm <- NA
    }
    
    else{
      hl_lm <- -1/coef(fit_lm)[[2]]
    }
    
    

    results <- rbind(results, c(hl_nls, hl_lm))
  }#end for loop
  names(results) <- c("hl_nls","hl_lm")
  results <- cbind(DGEList$genes, results)
  return(results)
}#end function