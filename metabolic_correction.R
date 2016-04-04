#' Normalizes metabolite data
#' 
#' Normalizes the intensity values of metabolite samples using the Virtual-QC method for systematic error correction
#'
#' @param samples a vector of the samples to be corrected
#' @param qcs a vector of the QCs to be used in correction
#' @param qc_positions a vector of the positions of each QC in the overall data;
#'        it is assumed that the last QC is the last value in the overall data.
#' @param w the number of adjacent QCs to use in creating the regression line for each sample in the Virtual-QC method
#'
#' @return A vector of normalized intensity values
#' @export
#'
#' @examples
normalize <- function(samples, qcs, qc_positions, w){
  n <- max(qc_positions)
  sample_positions <- setdiff(1:n, qc_positions)
  normalized <- mapply(sys_correct, samples, sample_positions, 
                       MoreArgs=list(qcs=qcs, qc_positions=qc_positions, w=w))
  return(normalized)
}

#' Systematically corrects a metabolite sample
#' 
#' Systematically corrects the intensity value of a metabolite sample using the Virtual-QC method
#'
#' @param sample the sample intensity value to be corrected
#' @param i the sample's position in the overall ordering, including QCs
#' @param qcs a vector of the QCs to be used in correction
#' @param qc_positions a vector of the positions of each QC in the overall data
#' @param w the number of adjacent QCs to use in creating the regression line for each sample in the Virtual-QC method
#'
#' @return A normalized intensity value
#' @export
#'
#' @examples
sys_correct <- function(sample, i, qcs, qc_positions, w){
  #calculate the index of the QC this sample comes after
  after_qc <- 1
  while(qc_positions[after_qc+1] < i) after_qc = after_qc + 1
  from <- max(1, after_qc - (w %/% 2) + 1)
  to <- min(after_qc + (w %/% 2), length(qcs))
  cf <- get_correction_factor(qcs, qc_positions, from, to, i)
  corrected <- sample / cf
  return(corrected)
}

#' Calculates a correction factor
#' 
#' Calculates a correction factor to use for correcting an intensity value in the Virtual-QC method
#'
#' @param qcs a vector of the QCs to be used in correction
#' @param qc_positions a vector of the positions of each QC in the overall data
#' @param from the smallest QC index to be considered in the regression line
#' @param to the largest QC index to be considered in the regression line
#' @param i the sample's position in the overall ordering
#'
#' @return The correction factor for sample i
#' @export
#'
#' @examples
get_correction_factor <- function(qcs, qc_positions, from, to, i){
  regression <- lm(qcs[from:to] ~ qc_positions[from:to])
  a <- regression$coefficients[[2]]
  b <- regression$coefficients[[1]]
  return(a*(i+1) + b)
}

#' Corrects QCs using regression lines
#' 
#' Corrects QC intensity values by creating regression lines of w adjacent QCs 
#'   and adjusting them to lie on those lines
#'
#' @param qcs the vector of QCs to be corrected
#' @param qc_positions a vector of the positions of each QC in the overall data
#' @param w the number of adjacent QCs to use in the regression line
#'
#' @return A vector of corrected QC intensity values
#' @export
#'
#' @examples
correct_qcs <- function(qcs, qc_positions, w){
  corrected <- sapply(1:length(qcs), correct_qc, qcs, qc_positions, w)
  return(corrected)
}

#' Corrects a QC using a regression line
#' 
#' Corrects a QC intensity value by adjusting it to lie on a regression line
#'   of w adjacent QCs
#'
#' @param i the index of the QC to be corrected
#' @param qcs the vector of QCs to be corrected
#' @param qc_positions a vector of the positions of each QC in the overall data
#' @param w the number of adjacent QCs to use in the regression line
#'
#' @return The corrected intensity value for QC i
#' @export
#'
#' @examples
correct_qc <- function(i, qcs, qc_positions, w){
  
  #TODO: this logic is stupid, make it better
  # index of QCs in overall data
  if(w %% 2 != 0){
    from <- max(0, i-((w-1)/2))
    to <- min(i+((w-1)/2), length(qcs))
    # if (i > w/2){
      # from_i <- i-((w-1)/2)
      # to_i <- min(i+((w-1)/2), length(qc_positions))
      
    # }
    # else {
      # from_i <- 1
      # to_i <- min(i+((w-1)/2), length(qc_positions))
    # }
    # from_i <- max(1, i-((w-1)/2))
    # to_i <- min(i+((w-1)/2), length(qc_positions))
  }
  else {
    if (i > (w/2)){
      from <- i-(w/2)
      to <- min(i+(w/2)-1, length(qcs))
      # from_i <- i-(w/2)
      # to_i <- min(i+(w/2)-1, length(qc_positions))
    }
    else {
      from <- 0
      to <- min(i+(w/2), length(qcs))
      # from_i <- 1
      # to_i <- min(i+(w/2), length(qc_positions))
    }
  }
  # from <- qc_positions[from_i]
  # to <- qc_positions[to_i]
  # print(paste("i=",i))
  # print(paste("i from", from_i, "to", to_i))
  # print(paste("x from", from, "to", to))
  
  # fit <- lm(unlist(combined_data[from:to])~seq(from, to))
  fit <- lm(qcs[from:to]~qc_positions[from:to])
  a <- fit$coefficients[[2]]
  b <- fit$coefficients[[1]]

  corrected <- a*(qc_positions[i]) + b
  # plot(qc_positions, qcs)
  # abline(fit)
  # points(qc_positions[i], corrected, col="blue", pch=19)
  # print(paste("x=",qc_positions[i]))
  # print(fit)
  # print(data[from:to])
  # print(seq(from, to))
  # plot(data[from:to])
  # abline(fit)
  # print(c(qc_positions[i], corrected))
  # points(list(x=positions[i], y=corrected), col="red", pch=19)

  return(corrected)
}

#' Corrects gross error in QCs
#' 
#' Corrects gross error in QC intensity values by multiplying by the ratios of two adjacent QCs
#'
#' @param qcs the vector of QCs to be corrected
#'
#' @return a vector of corrected QC intensity values
#' @export
#'
#' @examples
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

#' Calculates relative standard deviation
#' 
#' Calculates the relative standard deviation (standard deviation/mean) of a vector of values
#'
#' @param nums a vector of values
#'
#' @return The relative standard deviation of the values
#' @export
#'
#' @examples
relative_standard_deviation<- function(nums){
  n_mean = mean(nums)
  stand_dev <- sd(nums)
  rsd <- stand_dev / n_mean
  return(rsd)
}

#' Calculates ratios between QCs
#' 
#' Calculates the ratios between each adjacent pair of QC values, i.e. (QCi/QCi+1, QCi+1/QCi+2...)
#'
#' @param qcs a vector of QCs of size n
#'
#' @return A vector of ratios of size n-1
#' @export
#'
#' @examples
get_ratios <- function(qcs){
  ratios <- sapply(1:(length(qcs)-1), function(i){
    return(qcs[i]/qcs[i+1])
  })
  return(ratios)
}
