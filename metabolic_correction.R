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
  # look at plots of cf / datum, possibly using i
  return(corrected)
}

normalize <- function(data, qcs, betweens, qcpositions, w){
  indices = seq(1, length(data))
  # print(paste(length(data), length(betweens), length(indices), length(qcs)))
  normalized <- mapply(sys_correct, data, betweens, indices, MoreArgs=list(qcs=qcs, qcpositions=qcpositions, w=w))
  return(normalized)
}

correct_qcs <- function(qcs, positions, w){
  corrected = sapply(1:length(qcs), correct_qc, qcs, positions, w)
  # print(corrected)
  # return(listToPoints(corrected))
  return(corrected)
}

correct_qc <- function(i, qcs, positions, w){
  
  #TODO: this logic is stupid, make it better
  if(w %% 2 != 0){
    from <- max(0, i-((w-1)/2))
    to <- min(i+((w-1)/2), length(qcs))
  }
  else {
    if (i > (w/2)){
      from <- i-(w/2)
      to <- min(i+(w/2)-1, length(qcs))
    }
    else {
      from <- 0
      to <- min(i+(w/2), length(qcs))
    }
  }
  fit <- lm(qcs[from:to]~positions[from:to])
  a <- fit$coefficients[[2]]
  b <- fit$coefficients[[1]]
  
  corrected <- a*(positions[i]) + b
  # print(paste("x=",positions[i]))
  # print(fit)
  # plot(positions, qcs)
  # abline(fit)
  # print(c(positions[i], corrected))
  # points(list(x=positions[i], y=corrected), col="red", pch=19)

  return(corrected)
}

listToPoints <- function(l){
  p <- list(x=list(), y=list())
  for (point in l) {
    c(p$x, point$x)
    c(p$y, point$y)
  }
    
  return(p)
}


