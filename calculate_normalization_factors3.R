#Use non transformed data set fit nls

library(nls2)

calculate_normalization_factors3 <- function(DGEList = obj, method = c("average", "median", "peak")){
  
  #get the data, only keep the complete rows for calculation of correction coefficients
  obj <- DGEList
  data <- DGEList$cpm_trimmed #uncomment this line and comment next to inverse the order of trimming/ normalization
  #data <- cpm(obj)
  data <- data[complete.cases(data),]

  ##create a time point vector that fits the main data
  t <- sapply(strsplit(colnames(data),"t"), function(x){x[!x ==""]})
  
  ##loop over the genes to use for normalization and collect normalization factors
  ##normalization factors are coefficients that transform the curve into a perfect logarithmic decay
  corr_coeffs <- as.data.frame(data[0,])
  slopes <- numeric()
  for(i in 1:nrow(data)){
    # adjust the length of the time point vector to the actual data
    names(corr_coeffs) <- t
    t <- as.numeric(t)
    realt <- suppressWarnings(t+data[i,]-data[i,])
    
    ##Normalize the data so that t0 = 1
    values <- data[i,]/rep(data[i,1],length(realt))
    #values <- data[i,]
    #print(values)
    
    #fit
    #fit = tryCatch(nls(values ~ a*exp(-b*realt), data=list(realt,values), start=list(a= data[i,1],b=0.1)), error=function(e){print(paste("fitting failed row:" ,i))})
    #fit = tryCatch(nls(values ~ a*exp(-b*realt), data=list(realt,values), start=list(a= 1,b=0.1)), error=function(e){print(paste("fitting failed row:" ,i))})

    
    #Try fitting using nls2 (different starting values + optimization)
    
    st1 <- expand.grid(b = seq(0.05, 0.3, len = 10), a = seq(1, 1, len = 1))
    fit = tryCatch(nls2(values ~ a*exp(-b*realt), data=list(realt,values), start=st1, algorithm = "brute-force"), error=function(e){print(paste("fitting failed row:" ,i))})
    #fit = nls(values ~ a*exp(-b*realt), start = coef(fit))
    
    
    ifelse(grepl("fitting failed",fit), slope <- NA, slope <- coef(fit)[["b"]])
    #print(coef(fit))
    #exclude correction coefficients from negative slopes or failed fitting
    if(grepl("fitting failed",fit)){
      corr_coeffs <- rbind(corr_coeffs, values*NA)
      slopes <- c(slopes,NA)

    }
    
    else if(slope < 0){
      slopes <- c(slopes,NA)
    }
    else{
      corr_coeffs <- rbind(corr_coeffs,fitted(fit)/values)
      slopes <- c(slopes, slope)
    }
    
   
    }
    

 
  
  # write new normalization factors in the object
  corr_coeffs_no_NA <- corr_coeffs[complete.cases(corr_coeffs),]
  if(method == "mean"){
    DGEList$samples$norm.factors <- apply(corr_coeffs_no_NA, 2, mean)
  }
  if(method == "median"){
    DGEList$samples$norm.factors <- apply(corr_coeffs_no_NA, 2, median)
  }
  if(method == "peak"){
    DGEList$samples$norm.factors <- apply(corr_coeffs_no_NA, 2, function(z){
      den <- density(z)
      den$x[which.max(den$y)]
    } #closes function(z)
    )#close apply
  }#closes if
  
  print(paste("The calculated normalization factors were : "))
  print(DGEList$samples$norm.factors)
  
  ###PLots the distribution of the normalization factors for each time point
  par(mfrow = c(3,2), main = method,oma = c(0, 0, 2, 0))
  for(i in 1:length(DGEList$samples$norm.factors)){
    plot(density(corr_coeffs_no_NA[,i]), main = paste(" t =", colnames(corr_coeffs_no_NA)[i]))
  }
  mtext(paste("Distribution of correction coefficients"), outer = TRUE, cex = 1.5)

  DGEList$cpm_trimmed_normalized <- t(t(DGEList$cpm_trimmed)*DGEList$samples$norm.factors)
  
  return(DGEList)
  
  


}#closes function









