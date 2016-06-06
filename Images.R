make_rsd_plots_of_qcs <- function(d, start, positions, title){
  original_qcs <- get_all_qcs(d, start)
  
  plot(density(sapply(1:length(original_qcs), function(i){
    return(relative_standard_deviation(original_qcs[[i]]))
  } )), main=paste("RSDs of Original QCs", title))
  
  corrected_qcs <- lapply(1:length(corrected_qcs), function(i) {
    return(ratio_gross_correct(corrected_qcs[[i]]))
  })
  
  print(get_rsds_under_02(original_qcs))
  print(get_rsds_under_02(corrected_qcs))

  corrected_qcs <- lapply(corrected_qcs, regression_gross_correct,
                          w=5, qc_positions=positions)
  
  print(get_rsds_under_02(corrected_qcs))
  
  plot(density(sapply(1:length(corrected_qcs), function(i){
    return(relative_standard_deviation(corrected_qcs[[i]]))
  } )), main=paste("RSDS of Corrected QCs", title))
}

get_rsds_under_02 <- function(qcs){
  rsds <- sapply(1:length(qcs), function(i) {
    return(relative_standard_deviation(qcs[[i]]))
  })
  print(paste("max", max(rsds)))
  return(sum(rsds <= .2))
}

make_pca_images <- function(pca.data){
  pca <- prcomp(data, center=TRUE, scale.=TRUE)
  autoplot(pca, data=pca.data, colour='Type')
}