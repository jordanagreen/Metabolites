reg <- function(qcs, from, to){
  # print(paste("from ", 1, " to ", to-from+1))
  # print(paste(from, " ", to - from + 1))
  print(paste("from ", from, " to ", to))
  return(lm(qcs[from:to] ~ seq(1, 2)))
}

get_correction_factor <- function(qcs, from, to, i){
  #from <- max(0, i-w)
  #to <- min(i+w, length(qcs))
  #print(paste("from ", from, " to ", to))
  # regression = reg(qcs, from, to)
  # print(paste(from, to))
  # print(paste(qcs[from], qcs[to]))
  regression <- lm(qcs[from:to] ~ seq(1, to-from+1))
  # print(regression)
  a <- regression$coefficients[[2]]
  b <- regression$coefficients[[1]]
  #cat("a", a)
  #cat("b", b)
  return(a*(i+1) + b)
}

correct <- function(datum, between, i, qcs){
  # print(paste(datum, between, i))
  # print(qcs)
  cf <- get_correction_factor(qcs, between, between+1, i)
  corrected <- datum / cf
  # print(corrected)
  return(corrected)
}

normalize <- function(data, qcs, betweens){
  # scale <- 1:length(qcs)
  # normalized <- vector(,length(data))
  #print(scale)
  # for (i in 1:length(data)) {
  #   cf <- get_correction_factor(qcs, scale, between[i], between[i]+1)
  #   corrected <- data[i] / cf
  #   normalized[i] <- corrected
  # }
  indices = seq(1, length(data))
  normalized <- mapply(correct, data, betweens, indices, MoreArgs=list(qcs=qcs))
  return(normalized)
}

#add all data using qcs
#use mapply instead of for loop