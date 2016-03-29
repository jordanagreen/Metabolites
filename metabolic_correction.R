get_correction_factor <- function(qcs, qcpositions, from, to, i){
  regression <- lm(qcs[from:to] ~ qcpositions[from:to])
  a <- regression$coefficients[[2]]
  b <- regression$coefficients[[1]]
  return(a*(i+1) + b)
}

# get_correction_factor <- function(data, qcpositions, qc_i, w){
#   # index of qcs for window
#   from_i <- max(1, qc_i - (w %/% 2) + 1)
#   to_i <- min(qc_i + (w %/% 2), length(data))
#   from <- qcpositions[from_i]
#   to <- qcpositions[to_i]
#   regression <- lm(unlist(data[from:to])~seq(from, to))
#   a <- regression$coefficients[[2]]
#   b <- regression$coefficients[[1]]
#   return(a*(qc_i) + b)
# }

sys_correct <- function(sample, between, i, qcs, qcpositions, sample_positions, w){
  from <- max(1, between - (w %/% 2) + 1)
  to <- min(between + (w %/% 2), length(qcs))
  cf <- get_correction_factor(qcs, qcpositions, from, to, i)
  corrected <- sample / cf
  return(corrected)
}

normalize <- function(samples, qcs, betweens, qcpositions, sample_positions, w){
  normalized <- mapply(sys_correct, samples, betweens, sample_positions, 
                       MoreArgs=list(qcs=qcs, qcpositions=qcpositions, sample_positions, w=w))
  return(normalized)
}

correct_qcs <- function(data, qcpositions, w){
  corrected = sapply(1:length(qcpositions), correct_qc, data, qcpositions, w)
  # print(corrected)
  return(corrected)
}

correct_qc <- function(i, data, qcpositions, w){
  # data is samples and qcs, all in order of insertion
  # qcpositions is overall position of qc in data
  # i is qc #
  # w is window size
  
  #TODO: this logic is stupid, make it better
  # index of QCs in overall data
  if(w %% 2 != 0){
    if (i > w/2){
      from_i <- i-((w-1)/2)
      to_i <- min(i+((w-1)/2), length(qcpositions))
    }
    else {
      from_i <- 1
      to_i <- min(i+((w-1)/2), length(qcpositions))
    }
    # from_i <- max(1, i-((w-1)/2))
    # to_i <- min(i+((w-1)/2), length(qcpositions))
  }
  else {
    if (i > (w/2)){
      from_i <- i-(w/2)
      to_i <- min(i+(w/2)-1, length(qcpositions))
    }
    else {
      from_i <- 1
      to_i <- min(i+(w/2), length(qcpositions))
    }
  }
  from <- qcpositions[from_i]
  to <- qcpositions[to_i]
  # print(paste("i=",i))
  # print(paste("i from", from_i, "to", to_i))
  # print(paste("x from", from, "to", to))
  
  fit <- lm(unlist(data[from:to])~seq(from, to))
  a <- fit$coefficients[[2]]
  b <- fit$coefficients[[1]]
  
  corrected <- a*(qcpositions[i]) + b
  # print(paste("x=",qcpositions[i]))
  # print(fit)
  # print(data[from:to])
  # print(seq(from, to))
  # plot(data[from:to])
  # abline(fit)
  # print(c(qcpositions[i], corrected))
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


