ExperimentalVsPredicted <- function(submission.files, exp.val, column.P){
  path.prediction <- submission.files
  column.P1=3
  column.P2=5
  ifelse(column.P == column.P1, column.P.real <- 3, column.P.real <- 4)
  
  blueColor="#0000FFFF"

  #get the real data
  real.data <- exp.val
  realValue <- real.data[1:nrow(real.data), column.P.real]
  # axis dimension is defined by min aand max value of all predictions
  ifelse(column.P == columnP1, file.name <- "MO", file.name <- "MO_WT")
  png(paste("./results/Figure_S2_ExperimentalVsPredicted_", file.name, ".png", sep = ''), width = 9, height = 9, units = 'in', res = 300)
  
  par(mar=c(0,5,3,0)+.1, mfrow = c(2, 3), pty = "s", srt = 45) # decide number of plots in each row of the picture and total number of lines
  # check min and max value of predicted to decide plot axis length
  min.vector=rep(1.0, length(path.prediction)) # create a vector for min of each prediction file
  max.vector=rep(1.0, length(path.prediction)) # create a vector for max of each prediction file
  for(i in 1:length(path.prediction)) { 
    CurrentSubmission <- sort(path.prediction)[i]
    prediction.table <- read.table(file=CurrentSubmission, sep="\t", header=TRUE)
    predicted.value <- prediction.table[1:nrow(prediction.table), column.P]
    min.vector[i] <- min(predicted.value)
    max.vector[i] <- max(predicted.value)
  }
  min.max.predictions <- c(min(min.vector)*0.90, max(max.vector)*1.10) # vector for plot size
  
  #load the prediction
  for(i in 1:length(path.prediction)) { 
    CurrentSubmission <- sort(path.prediction)[i]
    prediction.table <- read.table(file=CurrentSubmission, sep="\t", header=TRUE)
    predicted.value <- prediction.table[1:nrow(prediction.table), column.P]
    submission.name <- gsub("-prediction_file-", ".", gsub("^Group_", "Submission ", gsub(".txt", "", basename(CurrentSubmission))))
    plot(min.max.predictions, min.max.predictions, main=submission.name, type= 'n', col = blueColor, xlab="Experimental values", ylab="Predicted values", cex.axis=1.0, cex.lab=1.1, cex.main=1.3, cex = 2)
    points(realValue, predicted.value, pch=20, cex = 2, col = blueColor)
    text(realValue, predicted.value, labels = real.data[1:nrow(real.data), 2], col = blueColor, pos = 4, cex=1.3) # comment this line to remove point lables
    segments(min.max.predictions[1], min.max.predictions[1], min.max.predictions[2], min.max.predictions[2], col = blueColor, lty= 2, lwd = 2)
  }
  
  dev.off() 
}
