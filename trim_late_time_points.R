trim_late_time_points <- function(DGeList = obj, CPMcutoff = 100000){
  #replaces late time points with cpm lower than threshold with NA
  #replaces all time points later than a NA time point with NA
  #if there is a cpm increase of more than 20 % between two time points, turn into NA
  #stores the "trimmed cpm" data in obj$trimmed
  q <- ifelse(cpm(DGeList) < CPMcutoff, NA, cpm(DGeList))
  difference <- t(apply(q, 1, diff))
  diffproc <- difference/q[,-ncol(q)]

  q[,-1] <- ifelse(diffproc > 0.2, NA, q[,-1])
  q[t(apply(is.na(q), 1, cumsum)) > 0 ] <- NA

  DGeList$cpm_trimmed <- q
  DGeList
 }

