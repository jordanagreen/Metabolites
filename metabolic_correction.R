#' Usage:
#' 
#' Step 1: Gather Data
#' Read the metabolite data into some dataset such that for each one, you have 
#' the sample values, the QC values, and the positions of each QC within the 
#' order of injection.
#' 
#' For example, if a line has thirteen values: samples 1-10 and QCs A, B, and C in positions
#' 4, 8, and 13, the samples vector would be c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#' the qcs vector would be c(A, B, C), and the qc_positions vector would be 
#' c(4, 8, 18). How these are calculated depends on how the original data is formatted.
#' 
#' Calculate the cutoff values beyond which any ratios between adjacent QCs are considered
#' to be outliers. This can be done by observing a density plot of all the ratios or some 
#' other method.
#' 
#' Step 2: QC Gross Error Correction
#' Apply the gross_correct_qcs function to each row of QCs in the dataset, with the
#' calculated cutoff values.
#' 
#' Step 3: QC Systematic Error Correction
#' Apply the systematic_correct_qcs function to each row of corrected QCs with its 
#' corresponding row of QC positions.
#' 
#' Step 4: Sample Normalization
#' Apply the normlaize_samples function to each row of samples with its corresponding
#' row of corrected QCs and QC positions.



#' Corrects gross error in QCs
#' 
#' Corrects gross error in QC intensity values by multiplying by the ratios of 
#' two adjacent QCs. Values which produce rations that fall outside the range of
#' 1 +- cutoff will be considered extreme outliers and replaced by the average 
#' of their neighbors
#' 
#' @param qcs the vector of QCs to be corrected
#' @param cutoff_min the value which any QCs lower than it will be considered
#'   outliers
#' @param cutoff_max the value which any QCs greater than it will be considered
#'   outliers
#'   
#' @return a vector of corrected QC intensity values
#' @export
#' 
#' @examples
gross_correct_qcs <- function(qcs, cutoff_min=.75, cutoff_max=1.25){
  
  corrected <- qcs
  ratios <- get_ratios(corrected)
  ratios_within_five_percent <- get_ratios_within_five_percent(ratios)
  
  # do correction until the ratios aren't getting any better
  t <- 1
  repeat{
    t <- t+1
    # correct QCs 2 to n-1, going forward
    for (i in 1:(length(corrected)-2)){
      # if the ratio is past the cutoff, QCi+1 is an outlier and should be replaced
      # by the average of QCi and QCi+2
      if (ratios[i] < cutoff_min | ratios[i] > cutoff_max){
        corrected[i+1] <- mean(c(corrected[i], corrected[i+2]))
      }
      else {
        # calculate the average of ratios R_i_i+1 and R_i+1_i+2
        avg_correction <- mean(c(ratios[i], ratios[i+1]))
        avg_corrected_i <- corrected[i+1] * avg_correction
        avg_corrected_i1 <- corrected[i+2] * avg_correction
        # if the ratio of the corrected QCs is better, replace the old QCs with them
        avg_corrected_ratio <- avg_corrected_i / avg_corrected_i1
        if (abs(1 - avg_corrected_ratio) < abs(1 - ratios[i])){
          corrected[i] <- avg_corrected_i
          corrected[i+1] <- avg_corrected_i1
        }
      }
    }
    
    # correct QCn
    
    # check if the last ratio is past the cutoff, and if so replace the last QCn
    # with the average of QCn-1 and QCn-2
    if (ratios[length(ratios)] < cutoff_min | ratios[length(ratios)] > cutoff_max){
      corrected[length(corrected)] <- mean(c(corrected[length(corrected)-1],
                                             corrected[length(corrected)-2]))
    }
    
    # correct QCs n-1 to 2, going backward
    for (i in length(corrected):3){
      if (ratios[i-1] < cutoff_min | ratios[i-1] > cutoff_max){
        corrected[i-1] <- mean(c(corrected[i], corrected[i=2]))
      }
      else {
        avg_correction <- mean(c(ratios[i-1], ratios[i-2]))
        avg_corrected_i <- corrected[i-1] * avg_correction
        avg_corrected_i1 <- corrected[i-2] * avg_correction
        avg_corrected_ratio <- avg_corrected_i / avg_corrected_i1
        if (abs(1 - avg_corrected_ratio) < abs(1 - ratios[i-1])){
          corrected[i] <- avg_corrected_i
          corrected[i-1] <- avg_corrected_i1
        }
      }
    }
    
    # correct QC1
    if (ratios[1] < cutoff_min | ratios[1] > cutoff_max){
      corrected[1] <- mean(c(corrected[2], corrected[3]))
    }
    # }
    # if we're still seeing improvement, keep going
    ratios <- get_ratios(corrected)
    new_ratios_within_five_percent <- get_ratios_within_five_percent(ratios)
    if (new_ratios_within_five_percent > ratios_within_five_percent){
      ratios_within_five_percent <- new_ratios_within_five_percent
      next
    }
    break
  }
  
  return(corrected)
}

#' Corrects systematic error in QCs using regression lines
#' 
#' Corrects QC intensity values by creating regression lines of w adjacent QCs 
#' and adjusting them to lie on those lines
#' 
#' @param qcs the vector of QCs to be corrected
#' @param qc_positions a vector of the positions of each QC in the overall data
#' @param w the number of adjacent QCs to use in the regression line
#'   
#' @return A vector of corrected QC intensity values
#' @export
#' 
#' @examples
systematic_correct_qcs <- function(qcs, qc_positions, w){
  corrected <- sapply(1:length(qcs), systematic_correct_qc, qcs, qc_positions, w)
  return(corrected)
}


#' Normalizes metabolite data
#' 
#' Normalizes the intensity values of metabolite samples using the Virtual-QC 
#' method for systematic error correction
#' 
#' @param samples a vector of the samples to be corrected
#' @param qcs a vector of the QCs to be used in correction
#' @param qc_positions a vector of the positions of each QC in the overall data;
#'   it is assumed that the last QC is the last value in the overall data.
#' @param w the number of adjacent QCs to use in creating the regression line 
#'   for each sample in the Virtual-QC method
#'   
#' @return A vector of normalized intensity values
#' @export
#' 
#' @examples
#' 
normalize_samples <- function(samples, qcs, qc_positions, w){
  n <- max(qc_positions)
  sample_positions <- setdiff(1:n, qc_positions)
  normalized <- mapply(systematic_correct_sample, samples, sample_positions, 
                       MoreArgs=list(qcs=qcs, qc_positions=qc_positions, w=w))
  return(normalized)
}

#' Systematically corrects a metabolite sample
#' 
#' Systematically corrects the intensity value of a metabolite sample using the
#' Virtual-QC method
#' 
#' @param sample the sample intensity value to be corrected
#' @param i the sample's position in the overall ordering, including QCs
#' @param qcs a vector of the QCs to be used in correction
#' @param qc_positions a vector of the positions of each QC in the overall data
#' @param w the number of adjacent QCs to use in creating the regression line
#'   for each sample in the Virtual-QC method
#'   
#' @return A normalized intensity value
#' @export
#' 
#' @examples
systematic_correct_sample <- function(sample, i, qcs, qc_positions, w){
  #calculate the index of the QC this sample comes after
  after_qc <- 1
  while(qc_positions[after_qc+1] < i) after_qc = after_qc + 1
  from <- max(1, after_qc - (w %/% 2) + 1)
  to <- min(after_qc + (w %/% 2), length(qcs))
  cf <- get_correction_factor(qcs, qc_positions, from, to, i)
  corrected <- (sample^2) / cf
  return(corrected)
}

#' Calculates a correction factor
#' 
#' Calculates a correction factor to use for correcting an intensity value in
#' the Virtual-QC method
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



#' Corrects a QC using a regression line
#' 
#' Corrects a QC intensity value by adjusting it to lie on a regression line of
#' w adjacent QCs
#' 
#' @param i the index of the QC to be corrected in the overall data
#' @param qcs the vector of QCs to be corrected
#' @param qc_positions a vector of the positions of each QC in the overall data
#' @param w the number of adjacent QCs to use in the regression line
#'   
#' @return The corrected intensity value for QC i
#' @export
#' 
#' @examples
systematic_correct_qc <- function(i, qcs, qc_positions, w){
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
  fit <- lm(qcs[from:to]~qc_positions[from:to])
  a <- fit$coefficients[[2]]
  b <- fit$coefficients[[1]]

  corrected <- a*(qc_positions[i]) + b
  return(corrected)
}



#' Calculates relative standard deviation
#' 
#' Calculates the relative standard deviation (standard deviation/mean) of a
#' vector of values
#' 
#' @param nums a vector of values
#'   
#' @return The relative standard deviation of the values
#' @export
#' 
#' @examples
relative_standard_deviation <- function(nums){
  n_mean = mean(nums)
  stand_dev <- sd(nums)
  rsd <- stand_dev / n_mean
  return(rsd)
}

#' Calculates ratios between QCs
#' 
#' Calculates the ratios between each adjacent pair of QC values, i.e.
#' (QCi/QCi+1, QCi+1/QCi+2...)
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

#' Calculates the amount of ratios between .95 and 1.05
#'
#' @param ratios a vector of ratios to check
#'
#' @return The total number of ratios between .95 and 1.05
#' @export
#'
#' @examples
get_ratios_within_five_percent <- function(ratios){
  within_five <- sum(ratios <= 1.05 & ratios >= .95)
  return(within_five)
}
