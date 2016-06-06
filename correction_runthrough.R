dataset <- read.csv("dataset.csv", header=F)
dataset2 <- read.csv("dataset-2.csv", header=F)
pca_dataset <- read.csv("pca-dataset.csv", header=T, row.names=1)
pca_dataset_2 <- read.csv("pca-dataset-2.csv", header=T, row.names=1)
cor_pca_dataset_2 <- read.csv("pca-dataset-2-corrected.csv", header=T, row.names=1)
qcstart <- 80
qcstart2 <- 81
qcpositions <- c(1, 12, 23, 34, 43, 44, 55, 66, 77, 88)
qcpositions2 <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 22, 33, 44, 55, 66, 77, 88, 98)

run <- function(d, start, qcpos){
  w <- 5
  n <- 20
  
  # original
  original_qcs <- lapply(1:nrow(d), get_qcs_by_metabolite, d=d, start=start)
  samples <- lapply(1:nrow(d), get_metabolite_data, d=d, start=start)
  # plot_rsds_and_save(original_qcs, "corrected-0a-original")
  
  corrected_qcs <- original_qcs

  # ratio gross correct
  for (i in 1:n){
    corrected_qcs <- lapply(1:length(corrected_qcs), function(i) {
      # print(corrected_qcs[[i]])
      return(ratio_gross_correct(corrected_qcs[[i]]))
    })
    # plot_rsds_and_save(corrected_qcs, paste("corrected", i, sep="-"))
    # print(paste("Done with correction", i))
  }
  return(corrected_qcs)
  
  # regression gross correct
  corrected_qcs <- lapply(corrected_qcs, regression_gross_correct,
                          w=w, qc_positions=qcpos)
  return(corrected_qcs)
  
  # plot_rsds_and_save(corrected_qcs, paste("corrected", n+1, sep=""))
  
  old_normalized <- mapply(normalize, samples, original_qcs, 
                           MoreArgs=list(qc_positions=qcpos, w=w), 
                           SIMPLIFY=FALSE)
  
  # normalize data with corrected QCs
  
  # normalized <- mapply(normalize, samples, corrected_qcs,
  #                      MoreArgs=list(qc_positions=qc_positions, w=w), 
  #                      SIMPLIFY=FALSE)
  print("Done normalizing")
  # return(normalized)
  # for (i in 1:length(normalized)){
  #   plot_metabolite_data(normalized[[i]], old_normalized[[i]], sample_positions, i)
  # }
  
  
  # corrected_rsds <- sapply(corrected_qcs, relative_standard_deviation)
  # normalized_rsds <- mapply(get_normalized_rsd, samples, corrected_qcs, 
                            # MoreArgs = list(w=w, qc_positions=qc_positions))
  # barplot(rsds)
  # barplot(corrected_rsds, main=paste("QCs Corrected", n, "times"))

  # original_normalized_rsds <- mapply(get_normalized_rsd, samples, original_qcs,
                                     # MoreArgs = list(w=w, qc_positions=qc_positions))
  # barplot(original_normalized_rsds, main=paste("Original normalized RSDs"))
  # barplot(normalized_rsds, main=paste("Normalized corrected", n, "times"))
  # print(paste(sum(original_normalized_rsds<=.2),"->", 
              # sum(normalized_rsds<=.2),"/", length(normalized_rsds)))
  # return(normalized_rsds)
  # return(sapply(corrected_qcs, relative_standard_deviation))
  
  # original_rsds <- sapply(original_qcs, relative_standard_deviation)
  # corrected_rsds <- sapply(corrected_qcs, relative_standard_deviation)
  # fixed_rsds <- mapply(function(o, c) return(o > c), original_rsds, corrected_rsds)
  # print(paste(sum(fixed_rsds), "/", length(fixed_rsds), "RSDs improved"))
  # 
  # return(lapply(corrected_qcs, get_ratios))
}

get_qcs_by_metabolite <- function(i, d, start){
  row <- d[i,]
  qcs <- as.numeric(row[start:length(row)])
  return(qcs)
}

get_all_qcs <- function(d, start){
  return(lapply(1:nrow(d), get_qcs_by_metabolite, d=d, start=start))
}

get_metabolite_data <- function(i, d, start){
  row <- d[i,]
  # print(row[2])
  data <- as.numeric(row[2:(start-1)])
  return(data)
}

get_all_metabolites <- function(d, start){
  return(t(as.data.frame(sapply((1:nrow(d)), get_metabolite_data, start=start, d=d))))
}

get_all_data_by_row <- function(d, i){
  row <- d[i,]
  return(as.numeric(row[2:length(row)]))
}

get_all_data <- function(d){
  all_data <- d
  return(all_data)
}

get_normalized_rsd <- function(samples, qcs, qc_positions, w){
  normalized <- normalize(samples, qcs, qc_positions, w)
  return(relative_standard_deviation(normalized))
}



plot_rsds_and_save <- function(qcs, f=""){
  rsds <- sapply(qcs, relative_standard_deviation)
  rsds_under_twenty_percent <- sum(rsds <= .2)
  ratios_within_five_percent <- sum(sapply(qcs, get_ratios_within_five_percent))
  png(filename=paste("img-progress/", f, ".png", sep=""))
  barplot(rsds, main=paste(rsds_under_twenty_percent, "RSDs <= .2,",
          ratios_within_five_percent, "ratios = 1+-.5"))
  abline(a=.2, b=0, col="red")
  dev.off()
}

get_all_ratios <- function(d, start){
  qcs <- lapply(1:nrow(d), get_qcs_by_metabolite, d=d, start=start)
  ratios <- unlist(lapply(qcs, get_ratios))
  return(ratios)
}

plot_metabolite_data <- function(samples, old_samples, sample_positions, i){
  png(filename=paste("img-progress/normalized/", i, ".png", sep=""))
  plot(sample_positions, samples, main=paste("Metabolite", i), col="blue", pch=19)
  points(sample_positions, old_samples, col="red", pch=19)
  dev.off()
}

correct_dataset <- function(d){
  normalized <- apply(d[-ncol(d)], 2, function(metabolite){
    samples <- metabolite[1:(qcstart-2)]
    qcs <- metabolite[(qcstart-1):length(metabolite)]
    # ratio gross correct
    for (i in 1:n){
      qcs <- ratio_gross_correct(qcs)
    }
    # regression gross correct
    qcs <- regression_gross_correct(qcs, qc_positions, 5)
    normalized <- normalize(samples, qcs, qc_positions, 3)
    return(c(normalized, qcs))
  })
  print("Done normalizing")
  # return(normalized)
  for (i in 1:(ncol(d)-1)){
    d[,i] <- normalized[,i]
  }
  return(d)
}

plot_correction <- function(i){
  qcs <- get_qcs_by_metabolite(i)
  plot(qcpositions, qcs, xlab="Order of Injection", ylab="Intensity")
  # title(main="Example of QC Correction")
  points(qcpositions, regression_gross_correct(ratio_gross_correct(qcs), qcpositions, 5), col="blue", pch=19)
  # legend("bottomleft", c("Original", "Corrected"), pch=c(1,19), col=c("black", "blue"))
}