get_correction_factor <- function(qcs, from, to, i){
  # print(class(qcs))
  regression <- lm(qcs[from:to] ~ seq(to, from))
  a <- regression$coefficients[[2]]
  b <- regression$coefficients[[1]]
  return(a*(i+1) + b)
}

sys_correct <- function(datum, between, i, qcs, w){
  from <- max(1, between - (w %/% 2) + 1)
  to <- min(between + 1 + (w %/% 2), length(qcs))
  cf <- get_correction_factor(qcs, from, to, i)
  corrected <- datum / cf
  return(corrected)
}

normalize <- function(data, qcs, betweens, w){
  indices = seq(1, length(data))
  # print(paste(length(data), length(betweens), length(indices), length(qcs)))
  normalized <- mapply(sys_correct, data, betweens, indices, MoreArgs=list(qcs=qcs, w=w))
  return(normalized)
}



