trim_late_time_points <- function(DGEList =obj, CPMcutoff= 0.5){
  #replaces late time points with cpm lower than threshold with NA
  #replaces all time points later than a NA time point with NA
  #if there is a cpm increase of more than 50 % between two time points, turn into NA
  #stores the "trimmed cpm" data in obj$trimmed
  obj <- DGEList
  q <- ifelse(cpm(obj)<CPMcutoff,NA,cpm(obj))
  q <- apply(q,1, function(x){
    for(i in 1:(length(x)-1)){
      if(((x[i+1]-x[i])/x[i]) > 0.5 & !is.na(x[i]) & !is.na(x[i+1])){
        x[i+1] <-NA
      }
      if(is.na(x[i])) {
        x[i+1] <- NA
      }
    }
    return(x)
  })

  obj$cpm_trimmed <- t(q)
  obj
  
  
  
  # q <- apply(q,1, function(x){
  #   for(i in 1:(length(x))){
  #     if(i < length(x)){
  #       if(((x[i+1]-x[i])/x[i]) > 0.2 & !is.na(x[i]) & !is.na(x[i+1])){
  #         x[i+1] <-NA
  #       }
  #     }
  #     if(is.na(x[i])) {
  #       x[i+1] <- NA
  #     }
  #   }
  #   return(x)
  # })
  # 
  # obj$cpm_trimmed <- t(q)
  # obj
  
}