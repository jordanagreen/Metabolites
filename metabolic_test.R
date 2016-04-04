betweens = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9)
dataset <- read.csv("dataset.csv", header=F)
qcstart <- 80
qc_positions <- c(1, 12, 23, 34, 43, 44, 55, 66, 77, 88)
sample_positions <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87)
combined_dataset <- read.csv("combined-dataset.csv", header=F)

get_data_and_qcs <- function(r){
  row <- combined_dataset[r,]
  return(row[2:length(row)])
}

testrow <- function(r){
  row <- dataset[r,]
  data <- as.numeric(row[3:qcstart-1])
  qcs <- as.numeric(row[qcstart:length(row)])
  par(mfrow = c(2, 1))
  plot(sample_positions, data, ylab="Intensity", xlab="Sample", main=paste("Original-", r, sep=""))
  plot(normalize(data, qcs, betweens, qc_positions, sample_positions, 2), main=paste("Normalized-", r, sep=""), ylab="Intensity", xlab="Sample")
}

test <- function(){
  data = c(1350290.99, 1383933.66, 1312871.32, 1310350.91, 1377561.05, 1372711.84, 1185227.79, 1293937.46, 1323975.87, 1140901.56, 1178429.53, 1338196.56, 1583288.71, 1449223.07, 1339895.09, 1311917.01, 1320321.76, 1310053.52, 1462071.76, 1400833.8, 1237759.06, 1291019.17, 958689.136, 1334708.01, 1189461.73, 1405293.91, 1062646.89, 1007465.82, 1394640.09, 1286606.02, 1222366.38, 1300354.64, 1133538.98, 1167052.53, 847428.638, 1254632.63, 980865.363, 1180603.83, 1179241.8, 1086125.52, 1315280.88, 1210457.53, 959763.963, 1094998.44, 1093956.32, 1181681.2, 1167093.43, 1006096.58, 1006781.36, 868498.665, 1196192.41, 1105860.67, 891654.625, 906416.079, 997145.024, 885163.246, 1053738.46, 928402.461, 812872.422, 848232.995, 853753.06, 800697.395, 1589040.16, 791363.152, 801715.379, 1141246.8, 818148.32, 807758.833, 934033.873, 839706.836, 874191.778, 1000508.67, 757952, 811035.543, 824572.941, 773340.875, 1743512.75, 725977.345)
  qcs = c(1441647.35, 1474580.65, 1403129.89, 1385726.05, 1201691.6, 1027367.38, 1231195.96, 1007911.46, 920991.688, 856119.891)
  betweens = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9)
  normalized = normalize(data, qcs, betweens, 2)
  print(normalized)
  plot(normalized)
  return(normalized)
}

plotAndSaveAll <- function(n){
  dataset <- read.csv("dataset.csv", header=F)
  lapply(1:nrow(dataset), plotAndSave, n)
}

plotAndSave <- function(r, n){
  row <- dataset[r,]
  qcstart <- 80
  
  data <- as.numeric(row[3:qcstart-1])
  # print(data)
  qcs <- as.numeric(row[qcstart:length(row)])
  png(filename=paste(n, "/", r, ".png", sep=""))
  par(mfrow = c(2, 1))
  # plot(data, ylab="Intensity", xlab="Sample", main=paste("Original-", r, sep=""))
  # points( c(0, 10, 20, 30, 30, 38, 48, 58, 68, 78), qcs, pch=21, bg="red")
  plot(qc_positions, qcs, pch=21, bg="red", main=paste("Original-", r, sep=""), ylab="Intensity", xlab="Sample")
  points(sample_positions, data)
  normalized <- normalize(data, qcs, betweens, qc_positions, sample_positions, 2)
  plot(sample_positions, normalized, main=paste("Normalized-", r, sep=""), ylab="Intensity", xlab="Sample")
  # plot(qcs, main=paste("QC-",r,sep=""), xlab="QC", ylab="Intensity")
  
  # print("done")
  dev.off()
}

only_plot <- function(r){
  par(mfrow = c(1, 1))
  row <- dataset[r,]
  # qcstart <- 80
  # betweens = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9)
  
  data <- as.numeric(row[3:qcstart-1])
  # print(data)
  qcs <- as.numeric(row[qcstart:length(row)])
  plot(sample_positions, data, ylab="Intensity", xlab="Sample", main=paste("Original-", r, sep=""))

  # plot(normalize(data, qcs, betweens, 2), main=paste("Normalized-", r, sep=""), ylab="Intensity", xlab="Sample")
  # plot(qcs, main=paste("QC-",r,sep=""), xlab="QC", ylab="Intensity")
  
  # print("done")
  # dev.off()
}

test_correct_qcs <- function(r){
  w <- 5
  par(mfrow = c(1, 1))
  
  # png(filename=paste("img-qcswithsamples/", r, ".png", sep=""))
  
  row <- dataset[r,]
  data <- as.numeric(row[3:qcstart-1])
  qcs <- as.numeric(row[qcstart:length(row)])
  corrected <- correct_qcs(qcs, qc_positions, 4)
  # c_row <- combined_dataset[r,]
  # combined <- as.numeric(c_row[2:length(c_row)])
  # corrected <- correct_qcs(combined, qc_positions, w)
  plot(sample_positions, data,  ylab="Intensity", xlab="Sample", main=paste("Metabolite", r))
  points(qc_positions, corrected, col="blue", pch=19)
  points(qc_positions, qcs, col="red", pch=19)
  abline(lm(qcs~qc_positions), col="red")
  abline(lm(corrected~qc_positions), col="blue")
  
  # dev.off()
}

test_plot_with_corrected_qcs <- function(r){
  w <- 5
  par(mfrow = c(1, 1))
  row <- dataset[r,]
  data <- as.numeric(row[3:qcstart-1])
  qcs <- as.numeric(row[qcstart:length(row)])
  # c_row <- combined_dataset[r,]
  # combined <- as.numeric(c_row[2:length(c_row)])
  # corrected <- correct_qcs(combined, qc_positions, w)
  corrected <- correct_qcs(qcs, qc_positions, w)
  normalized <- normalize(data, qcs, betweens, qc_positions, sample_positions, w)
  corrected_normalized <- normalize(data, corrected, betweens, qc_positions, sample_positions, w)
  plot(normalized, col="red", pch=19)
  points(corrected_normalized, col="blue", pch=19)
}


test_corrected_ratio <- function(r){
  row <- dataset[r,]
  data <- as.numeric(row[3:qcstart-1])
  qcs <- as.numeric(row[qcstart:length(row)])
  # c_row <- combined_dataset[r,]
  # combined <- as.numeric(c_row[2:length(c_row)])
  # corrected <- correct_qcs(combined, qc_positions, 3)
  corrected <- correct_qcs(qcs, qc_positions, w)
  ratios <- numeric(length(qcs)-1)
  for (i in 1:length(ratios)) {
    ratios[i] <- corrected[i]/corrected[i+1]
  }
  return(ratios)
}

test_ratios <- function(){
  ratios <- lapply(1:nrow(dataset), test_corrected_ratio)
  df <- data.frame(matrix(unlist(ratios), nrow=length(ratios), byrow=T))
  # results <- table(ratios)
  write.csv(df, file="ratios_with_samples.csv")
  return(df)
}

test_rsd <- function(r){
  row <- dataset[r,]
  qcs <- as.numeric(row[qcstart:length(row)])
  return(relative_standard_deviation(qcs))
}

graph_qc_correction <- function(r){
  row <- dataset[r,]
  data <- as.numeric(row[3:qcstart-1])
  qcs <- as.numeric(row[qcstart:length(row)])
  plot(qc_positions, qcs, col="black", pch=19)
  points(qc_positions, gross_correct(qcs), col="red", pch=19)
  points(qc_positions, correct_qcs(qcs, qc_positions, 5), col="green", pch=19)
  points(qc_positions, correct_qcs(gross_correct(qcs), qc_positions, 5), col="blue", pch=19)
  points(sample_positions, data)
}

graph_rsds <- function(){
  w <- 3
  # rsds <- apply(dataset, 1, function(x){return(relative_standard_deviation(x[2:length(x)]))})
  rsds <- sapply(1:nrow(dataset), function(r){
    row <- dataset[r,]
    data <- as.numeric(row[3:qcstart-1])
    qcs <- as.numeric(row[qcstart:length(row)])
    # c_row <- combined_dataset[r,]
    # combined <- as.numeric(c_row[2:length(c_row)])
    # corrected <- correct_qcs(combined, qc_positions, w)
    corrected <- correct_qcs(qcs, qc_positions, w)
    return(relative_standard_deviation(corrected))
    })
  # print(rsds)
  barplot(rsds, main="Relative Standard Deviations", xlab="Sample", ylab="RSD", ylim=c(0,2))
  abline(a=.2, b=0, col="red")
}

test_normalize_with_gross_correction <- function(r){
  w <- 3
  n <- 5
  par(mfrow = c(1, 1))
  row <- dataset[r,]
  data <- as.numeric(row[3:qcstart-1])
  qcs <- as.numeric(row[qcstart:length(row)])
  # c_row <- combined_dataset[r,]
  # combined <- as.numeric(c_row[2:length(c_row)])
  # corrected <- qcs
  corrected <- gross_correct(qcs)
  lapply(1:n, function(i) {
    corrected <<- gross_correct(corrected)
    # print(i)
    # print(corrected)
    })
  
  # combined_with_corrected_qcs <- sapply(1:length(combined), function(i){
  #   if (i %in% qc_positions){
  #     return(corrected[match(i, qc_positions)])
  #   }
  #   else return(data[i])
  # })
  corrected <- correct_qcs(qcs, qc_positions, 5)
  
  normalized <- normalize(data, qcs, betweens, qc_positions, sample_positions, w)
  corrected_normalized <- normalize(data, corrected, betweens, qc_positions, sample_positions, w)
  par(mfrow=c(1,1))
  nrsd <- relative_standard_deviation(normalized)
  cnrsd <- relative_standard_deviation(corrected_normalized)
  # print(paste(nrsd, cnrsd))
  # png(filename=paste("img-gross/", r, ".png", sep=""))
  # plot(normalized, col="red", pch=19, main=paste(r, "RSD", nrsd, "->", cnrsd))
  # points(corrected_normalized, col="blue", pch=19)
  # dev.off()
  return(cnrsd <= nrsd)
}

test_normalize_with_gross_correction_all <- function(){
  fixed <- lapply(1:nrow(dataset), function(r) test_normalize_with_gross_correction(r))
  print(paste(sum(unlist(fixed)), "/", length(fixed)))
}

test_gross_correction_rsds <- function(){
  n <- 10
  all_qcs <- lapply(1:nrow(dataset), function(r){
    row <- dataset[r,]
    qcs <- as.numeric(row[qcstart:length(row)])
    return(qcs)
  })
  lapply(1:n, function(i){
    all_qcs <<- lapply(all_qcs, function(qcs){
      return(gross_correct(qcs))
    })
    rsds <- sapply(all_qcs, function(qcs){
      return(relative_standard_deviation(qcs))
    })
    # print(rsds[3])
    # barplot(rsds)
    # barplot(rsds, main=paste("Correction", i, "-", "Total RSD < .2:", sum(rsds<.2)), xlab="Sample", ylab="RSD", ylim=c(0,2))
    # abline(a=.2, b=0, col="red")
  })
  print("done")
}

test_gross_correction_ratios <- function(){
  n <- 20
  ratios <- lapply(1:nrow(dataset), function(r){
    row <- dataset[r,]
    qcs <- as.numeric(row[qcstart:length(row)])
    return(get_ratios(qcs))
  })
  corrected_ratios <- lapply(1:nrow(dataset), function(r){
    row <- dataset[r,]
    data <- as.numeric(row[3:qcstart-1])
    qcs <- as.numeric(row[qcstart:length(row)])
    corrected_qcs <- qcs
    for(i in 1:n){
      corrected_qcs <- gross_correct(corrected_qcs)
    }
    corrected_qcs <- correct_qcs(corrected_qcs, qc_positions, 5)
    return(get_ratios(corrected_qcs))
  })
  fixed <- mapply(function(r, cr) return(abs(1-cr) <= abs(1-r)), ratios, corrected_ratios)
  print(paste(sum(unlist(fixed)), "/", length(fixed)))
  worse <- mapply(function(r, cr) return(abs(1-cr) > abs(1-r)), ratios, corrected_ratios)
  print(paste(sum(unlist(worse)), "/", length(worse)))
  
}
