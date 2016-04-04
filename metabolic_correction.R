normalize <- function(samples, qcs, betweens, qcpositions, sample_positions, w){
  normalized <- mapply(sys_correct, samples, betweens, sample_positions, 
                       MoreArgs=list(qcs=qcs, qcpositions=qcpositions, sample_positions, w=w))
  return(normalized)
}

sys_correct <- function(sample, between, i, qcs, qcpositions, sample_positions, w){
  from <- max(1, between - (w %/% 2) + 1)
  to <- min(between + (w %/% 2), length(qcs))
  cf <- get_correction_factor(qcs, qcpositions, from, to, i)
  corrected <- sample / cf
  return(corrected)
}

get_correction_factor <- function(qcs, qcpositions, from, to, i){
  regression <- lm(qcs[from:to] ~ qcpositions[from:to])
  a <- regression$coefficients[[2]]
  b <- regression$coefficients[[1]]
  return(a*(i+1) + b)
}

# correct_qcs <- function(combined_data, qcpositions, w){
#   corrected = sapply(1:length(qcpositions), correct_qc, combined_data, qcpositions, w)
#   # print(corrected)
#   return(corrected)
# }

correct_qcs <- function(qcs, qcpositions, w){
  corrected <- sapply(1:length(qcs), correct_qc, qcs, qcpositions, w)
  return(corrected)
}

# correct_qc <- function(i, combined_data, qcpositions, w){
correct_qc <- function(i, qcs, qcpositions, w){
  
  # qcpositions is overall position of qc in data
  # i is qc #
  # w is window size
  
  #TODO: this logic is stupid, make it better
  # index of QCs in overall data
  if(w %% 2 != 0){
    from <- max(0, i-((w-1)/2))
    to <- min(i+((w-1)/2), length(qcs))
    # if (i > w/2){
      # from_i <- i-((w-1)/2)
      # to_i <- min(i+((w-1)/2), length(qcpositions))
      
    # }
    # else {
      # from_i <- 1
      # to_i <- min(i+((w-1)/2), length(qcpositions))
    # }
    # from_i <- max(1, i-((w-1)/2))
    # to_i <- min(i+((w-1)/2), length(qcpositions))
  }
  else {
    if (i > (w/2)){
      from <- i-(w/2)
      to <- min(i+(w/2)-1, length(qcs))
      # from_i <- i-(w/2)
      # to_i <- min(i+(w/2)-1, length(qcpositions))
    }
    else {
      from <- 0
      to <- min(i+(w/2), length(qcs))
      # from_i <- 1
      # to_i <- min(i+(w/2), length(qcpositions))
    }
  }
  # from <- qcpositions[from_i]
  # to <- qcpositions[to_i]
  # print(paste("i=",i))
  # print(paste("i from", from_i, "to", to_i))
  # print(paste("x from", from, "to", to))
  
  # fit <- lm(unlist(combined_data[from:to])~seq(from, to))
  fit <- lm(qcs[from:to]~qcpositions[from:to])
  a <- fit$coefficients[[2]]
  b <- fit$coefficients[[1]]

  corrected <- a*(qcpositions[i]) + b
  # plot(qcpositions, qcs)
  # abline(fit)
  # points(qcpositions[i], corrected, col="blue", pch=19)
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

gross_correct <- function(qcs){
  rsd <- relative_standard_deviation(qcs)
  # print(paste("RSD", rsd))
  # Samples with an RSD greater than 20% need to be corrected further
  if (rsd > .2){
    corrected <- qcs
    ratios <- get_ratios(corrected)
    # print("ratios")
    # print(ratios)
    
    #For i = 1 to len(qcs)-2
    lapply(1:(length(corrected)-2), function(i){
      avg_correction <- mean(c(ratios[i], ratios[i+1]))
      ai_corrected <- corrected[i+1] * avg_correction
      ai1_corrected <- corrected[i+2] * avg_correction
      corrected_ratio <-ai_corrected / ai1_corrected
      # if this ratio is closer to 1 than ratios[i], keep it
      # print(paste(ratios[i], "to", corrected_ratio))
      if (abs(1-corrected_ratio) < abs(1-ratios[i])){
        # print(paste(i, "changed from", qcs[i], "to", ai_corrected))
        corrected[i] <<- ai_corrected #needs to modify array in outer scope, so <<-
        corrected[i+1] <<- ai1_corrected
        # print(paste(i, qcs[i]))
      }
    })
    
    # fix the last two qcs
    avg_correction <- mean(c(ratios[length(ratios)], ratios[length(ratios)-1]))
    an_corrected <- corrected[length(corrected)] * avg_correction
    an1_corrected <- corrected[length(corrected)-1] * avg_correction
    corrected_ratio <- an_corrected / an1_corrected
    if (abs(1-corrected_ratio) < abs(1-ratios[length(ratios)])){
      corrected[length(corrected)] <- an_corrected
      corrected[length(corrected)-1] <- an1_corrected
    }
    new_rsd <- relative_standard_deviation(corrected)
    if (new_rsd <= rsd){
      return(corrected)
    }
  }
  return(qcs)
}

relative_standard_deviation<- function(qcs){
  qc_mean = mean(qcs)
  stand_dev <- sd(qcs)
  rsd <- stand_dev / qc_mean
  return(rsd)
}

get_ratios <- function(qcs){
  ratios <- sapply(1:(length(qcs)-1), function(i){
    return(qcs[i]/qcs[i+1])
  })
  return(ratios)
}
