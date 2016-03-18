get_correction_factor <- function(qcs, qcpositions, from, to, i){
  # print(class(qcs))
  # print(paste(from, to))
  # print(paste(qcs[from:to], c(qcpositions[from], qcpositions[to])))
  regression <- lm(qcs[from:to] ~ c(qcpositions[from], qcpositions[to]))
  a <- regression$coefficients[[2]]
  b <- regression$coefficients[[1]]
  return(a*(i+1) + b)
}

sys_correct <- function(datum, between, i, qcs, qcpositions, w){
  from <- max(1, between - (w %/% 2) + 1)
  to <- min(between + (w %/% 2), length(qcs))
  # print(paste("between", between))
  # print(paste("from", from, "to", to))
  cf <- get_correction_factor(qcs, qcpositions, from, to, i)
  corrected <- datum / cf
  return(corrected)
}

normalize <- function(data, qcs, betweens, qcpositions, w){
  indices = seq(1, length(data))
  # print(paste(length(data), length(betweens), length(indices), length(qcs)))
  normalized <- mapply(sys_correct, data, betweens, indices, MoreArgs=list(qcs=qcs, qcpositions=qcpositions, w=w))
  return(normalized)
}

correct_qcs <- function(qcs, positions){
  fit <- lm(qcs~positions)
  a <- fit$coefficients[[2]]
  b <- fit$coefficients[[1]]
  corrected = mapply(qcs, positions, MoreArgs=list(correct_qc, a, b))
}

correct_qc <- function(qc, i, a, b){
  cf <- a*(i+1) + b
  corrected = qc*(1/cf)
  return(corrected)
}


