trim_late_time_points <- function(DGeList =obj, CPMcutoff= 0.5){
  #replaces late time points with cpm lower than threshold with NA
  #replaces all time points later than a NA time point with NA
  #if there is a cpm increase of more than 20 % between two time points, turn into NA
  #stores the "trimmed cpm" data in obj$trimmed
  q <- ifelse(cpm(obj)<CPMcutoff,NA,cpm(obj))
  q <- apply(q,1, function(x){
    for(i in 1:(length(x)-1)){
      if(((x[i+1]-x[i])/x[i]) > 0.2 & !is.na(x[i]) & !is.na(x[i+1])){
        x[i+1] <-NA
      }
      if(is.na(x[i])) {
        x[i+1] <- NA
      }
    }
    return(x)
  })
  
  obj$trimmed <- t(q)
  obj
  
  
}