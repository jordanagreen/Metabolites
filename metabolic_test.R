betweens = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9)
dataset <- read.csv("dataset.csv", header=F)
qcstart <- 80
qcpositions <- c(0, 10, 20, 30, 34, 38, 48, 58, 68, 78)

testrow <- function(r){
  
  # betweens = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9)
  row <- dataset[r,]
  data <- as.numeric(row[3:qcstart-1])
  qcs <- as.numeric(row[qcstart:length(row)])
  par(mfrow = c(2, 1))
  plot(data, ylab="Intensity", xlab="Sample", main=paste("Original-", r, sep=""))
  plot(normalize(data, qcs, betweens, qcpositions, 2), main=paste("Normalized-", r, sep=""), ylab="Intensity", xlab="Sample")
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
  plot(qcpositions, qcs, pch=21, bg="red", main=paste("Original-", r, sep=""), ylab="Intensity", xlab="Sample")
  points(data)
  normalized <- normalize(data, qcs, betweens, qcpositions, 2)
  plot(normalized, main=paste("Normalized-", r, sep=""), ylab="Intensity", xlab="Sample")
  # plot(qcs, main=paste("QC-",r,sep=""), xlab="QC", ylab="Intensity")
  
  # print("done")
  dev.off()
}

only_plot <- function(r){
  row <- dataset[r,]
  # qcstart <- 80
  betweens = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9)
  
  data <- as.numeric(row[3:qcstart-1])
  # print(data)
  qcs <- as.numeric(row[qcstart:length(row)])
  plot(data, ylab="Intensity", xlab="Sample", main=paste("Original-", r, sep=""))

  # plot(normalize(data, qcs, betweens, 2), main=paste("Normalized-", r, sep=""), ylab="Intensity", xlab="Sample")
  # plot(qcs, main=paste("QC-",r,sep=""), xlab="QC", ylab="Intensity")
  
  # print("done")
  # dev.off()
}

test_correct_qcs <- function(r){
  row <- dataset[r,]
  qcs <- as.numeric(row[qcstart:length(row)])
}