CorrelationHeatMap <- function(path.prediction)
{
  submission.name <- gsub("-prediction_file-", ".", gsub("^Group_", "Submission ", gsub(".txt$", "", basename(path.prediction))))
  correlation.p.value.MO = data.frame(matrix(NA, length(path.prediction), length(path.prediction)))
  rownames(correlation.p.value.MO) <- submission.name
  colnames(correlation.p.value.MO) <- submission.name
  correlation.p.value.MO.WT = data.frame(matrix(NA, length(path.prediction), length(path.prediction)))
  rownames(correlation.p.value.MO.WT) <- submission.name
  colnames(correlation.p.value.MO.WT) <- submission.name
  
  prediction.p.value.MO.array = data.frame(matrix(NA, 38, length(path.prediction)))
  prediction.p.value.MO.WT.array = data.frame(matrix(NA, 38, length(path.prediction)))
  colnames(prediction.p.value.MO.array) <- submission.name
  colnames(prediction.p.value.MO.WT.array) <- colnames(prediction.p.value.MO.array)
  
  for(i in 1:length(path.prediction)) {
    CurrentSubmission <- path.prediction[i]
    cat("****************** ", submission.name[i], " ******************\n")
    
    prediction.table <- read.table(file = CurrentSubmission, sep = "\t", quote = "", header = TRUE)
    prediction.p.value.MO.array[i] <- prediction.table$P.value_MO
    prediction.p.value.MO.WT.array[i] <- prediction.table$P.value_MO.WT
  }
  
  for(i in 1:length(path.prediction)){
    for(j in 1:length(path.prediction)){
      correlation.p.value.MO[i,j] <- cor(prediction.p.value.MO.array[i] , prediction.p.value.MO.array[j] , method = "pearson")
      correlation.p.value.MO.WT[i,j] <- cor(prediction.p.value.MO.WT.array[i] , prediction.p.value.MO.WT.array[j] , method = "pearson")
    }
  }
  
  
  file.name <- paste('./results/Figure_1_A_submissions_correlation_MO.png', sep="")
  png(file.name, width = 6, height = 6, units = 'in', res = 300)
  threshold = 0
  complete.score.matrix.sorted <- correlation.p.value.MO
  par(mar=c(0,2,5,0)+.1)
  cellcolors<-matrix(NA, ncol(complete.score.matrix.sorted), ncol(complete.score.matrix.sorted))
  cellcolors[!is.na(complete.score.matrix.sorted) & complete.score.matrix.sorted > threshold] <- color.scale(complete.score.matrix.sorted[complete.score.matrix.sorted > threshold & !is.na(complete.score.matrix.sorted)], cs1 = 0, cs2 = c(0.9, 0.4), cs3 = 0)
  cellcolors[!is.na(complete.score.matrix.sorted) & complete.score.matrix.sorted < threshold] <- color.scale(complete.score.matrix.sorted[complete.score.matrix.sorted < threshold & !is.na(complete.score.matrix.sorted)], cs1 = c(0.4, 0.9), cs2 = 0, cs3 = 0)
  cellcolors[!is.na(complete.score.matrix.sorted) & complete.score.matrix.sorted == threshold] <- "#FFFFFFFF"
  cellcolors[is.na(complete.score.matrix.sorted)] <- "#000000FF"
  plt <- color2D.matplot(complete.score.matrix.sorted, axes=FALSE, cellcolors=cellcolors, show.values = 2, xlab = "", ylab = "", main="Correlation among Submission MO")
  simple.name <- gsub("-prediction_file-", ".", gsub("^Group_", "", gsub(".txt$", "", basename(path.prediction))))
  labels <- simple.name
  axis(3, at=0.5:nrow(complete.score.matrix.sorted),  labels = labels)
  axis(2, at=0.5:nrow(complete.score.matrix.sorted), labels = rev(labels))
  dev.off()
  
  file.name <- paste('./results/Figure_1_B_submissions_correlation_MO_WT.png', sep="")
  png(file.name, width = 6, height = 6, units = 'in', res = 300)
  threshold = 0
  complete.score.matrix.sorted <- correlation.p.value.MO.WT
  par(mar=c(0,2,5,0)+.1)
  cellcolors<-matrix(NA, ncol(complete.score.matrix.sorted), ncol(complete.score.matrix.sorted))
  cellcolors[!is.na(complete.score.matrix.sorted) & complete.score.matrix.sorted > threshold] <- color.scale(complete.score.matrix.sorted[complete.score.matrix.sorted > threshold & !is.na(complete.score.matrix.sorted)], cs1 = 0, cs2 = c(0.9, 0.4), cs3 = 0)
  cellcolors[!is.na(complete.score.matrix.sorted) & complete.score.matrix.sorted < threshold] <- color.scale(complete.score.matrix.sorted[complete.score.matrix.sorted < threshold & !is.na(complete.score.matrix.sorted)], cs1 = c(0.4, 0.9), cs2 = 0, cs3 = 0)
  cellcolors[!is.na(complete.score.matrix.sorted) & complete.score.matrix.sorted == threshold] <- "#FFFFFFFF"
  cellcolors[is.na(complete.score.matrix.sorted)] <- "#000000FF"
  plt <- color2D.matplot(complete.score.matrix.sorted, axes=FALSE, cellcolors=cellcolors, show.values = 2, xlab = "", ylab = "", main="Correlation among Submission MO+WT")
  simple.name <- gsub("-prediction_file-", ".", gsub("^Group_", "", gsub(".txt$", "", basename(path.prediction))))
  labels <- simple.name
  axis(3, at=0.5:nrow(complete.score.matrix.sorted),  labels = labels)
  axis(2, at=0.5:nrow(complete.score.matrix.sorted), labels = rev(labels))
  dev.off()
}

ConfusionMatrixRandomSubmissions <- function(exp.val, path.prediction)
{
  Experimental.functional.effect <- exp.val$Functional.effect
  
  complete.score.matrix = data.frame(matrix(NA, length(path.prediction), 5))
  rownames(complete.score.matrix) <- paste("Random submission", seq(1:length(path.prediction)))
  colnames(complete.score.matrix) <- c("BACC", "MCC", "F1", "TPR", "TNR")
  
  for(i in 1:ncol(path.prediction)) {
    submission.name <- rownames(complete.score.matrix)[i]
    cat("****************** ", submission.name, " ******************\n")
    
    prediction.functional.effect <- path.prediction[,i]
    
    class = "0"
    TN.condition <- which(Experimental.functional.effect == class & prediction.functional.effect == class)
    FP.condition <- which(Experimental.functional.effect == class & prediction.functional.effect != class)
    FN.condition <- which(Experimental.functional.effect != class & prediction.functional.effect == class)
    TP.condition <- which(Experimental.functional.effect != class & prediction.functional.effect != class)
    TP <- length(TP.condition)
    FN <- length(FN.condition)
    FP <- length(FP.condition)
    TN <- length(TN.condition)
    
    submission.scores <- computeBasicStatistic(file.name=NULL, class=class, TP=TP, FN=FN, FP=FP, TN=TN)
    complete.score.matrix[i, ] <- submission.scores
  }
  
  return(complete.score.matrix)
}


ConfusionMatrix <- function(exp.val, path.prediction)
{
  Experimental.functional.effect <- exp.val$Functional.effect
  Experimental.label <- exp.val$Variant
  Experimental.p.value.from.MO <- exp.val$p.value.from.MO
  Experimental.p.value.from.MO.WT <- exp.val$p.value.from.MO.WT

  complete.score.matrix = data.frame(matrix(NA, length(path.prediction), 5))
  rownames(complete.score.matrix) <- gsub("-prediction_file-", ".", gsub("^Group_", "Submission ", gsub(".txt$", "", basename(path.prediction))))
  colnames(complete.score.matrix) <- c("BACC", "MCC", "F1", "TPR", "TNR")
  
  which.phenotyphe.matrix = array(0, dim=c(nrow(exp.val), 4, length(path.prediction)))
  rownames(which.phenotyphe.matrix) <- seq(1, nrow(exp.val), by=1)
  colnames(which.phenotyphe.matrix) <- c("TP", "FN", "FP", "TN")
  
  for(i in 1:length(path.prediction)) {
    CurrentSubmission <- sort(path.prediction)[i]
    submission.name <- gsub("-prediction_file-", ".", gsub("^Group_", "Submission ", gsub(".txt$", "", basename(CurrentSubmission))))
    cat("****************** ", submission.name, " ******************\n")
    
    prediction.table <- read.table(file = CurrentSubmission, sep = "\t", quote = "", header = TRUE)
    prediction.functional.effect <- prediction.table$Functional_effect
    prediction.p.value.MO <- prediction.table$P.value_MO
    prediction.SD.MO <- prediction.table$Standard_deviation
    prediction.p.value.MO.WT <- prediction.table$P.value_MO.WT
    prediction.SD.MO.WT <- prediction.table$Standard_deviation.1
    
    file.name <- paste('./results/confusion_matrix_balanced_', submission.name,'.txt', sep="")
    if(file.exists(file.name))
      file.remove(file.name)
    
    class = "0"
    TN.condition <- which(Experimental.functional.effect == class & prediction.functional.effect == class)
    FP.condition <- which(Experimental.functional.effect == class & prediction.functional.effect != class)
    FN.condition <- which(Experimental.functional.effect != class & prediction.functional.effect == class)
    TP.condition <- which(Experimental.functional.effect != class & prediction.functional.effect != class)
    TP <- length(TP.condition)
    FN <- length(FN.condition)
    FP <- length(FP.condition)
    TN <- length(TN.condition)
    
    submission.scores <- computeBasicStatistic(file.name=file.name, class=class, TP=TP, FN=FN, FP=FP, TN=TN)
    complete.score.matrix[i, ] <- submission.scores
    which.phenotyphe.matrix[TP.condition, 1, i] <- 1
    which.phenotyphe.matrix[FN.condition, 2, i] <- 1
    which.phenotyphe.matrix[FP.condition, 3, i] <- 1
    which.phenotyphe.matrix[TN.condition, 4, i] <- 1
  }
  
  computePValueSDPlotMultiRow(path.prediction=path.prediction, WT=FALSE, Experimental.p.value=Experimental.p.value.from.MO, complete=FALSE)
  computePValueSDPlotMultiRow(path.prediction=path.prediction, WT=TRUE, Experimental.p.value=Experimental.p.value.from.MO.WT, complete=FALSE)
  if(TRUE){
    computePValueSDPlotMultiRow(path.prediction=path.prediction, WT=FALSE, Experimental.p.value=Experimental.p.value.from.MO, complete=TRUE)
    computePValueSDPlotMultiRow(path.prediction=path.prediction, WT=TRUE, Experimental.p.value=Experimental.p.value.from.MO.WT, complete=TRUE)
  }
     
  write.table(round(complete.score.matrix, 2), file='./results/scores_indexes_by_submission.txt', sep = "\t" , quote = F, col.names = NA)
  
  return(list(complete.score.matrix, which.phenotyphe.matrix))
}

computeBasicStatistic <- function(file.name, class, TP, FN, FP, TN)
{
  PPV <- TP/(TP+FP)
  NPV <- TN/(TN+FN)
  TPR <- TP/(TP+FN)
  TNR <- TN/(TN+FP)
  MCC <- ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  ACC <- (TP+TN)/(TP+TN+FP+FN)
  BACC <- (TPR+TNR)/2
  F1 <- (2*PPV*TPR)/(PPV+TPR)
  confusion.matrix <- data.frame(c(TP, FN), c(FP, TN))
  rownames(confusion.matrix) <- c("Pred Benign", "Pred Others")
  colnames(confusion.matrix) <- c("Obs Benign", "Obs Others")
  print(confusion.matrix)
  cat(paste("\nPPV:", PPV, "\nNPV:", NPV, "\nTPR:", TPR, "\nTNR:", TNR, "\nMCC:", MCC, "\nACC:", ACC, "\nBACC:", BACC, "\nF1:", F1, "\n\n************************************\n"))
  if(!is.null(file.name)) 
  {
    write.table(confusion.matrix, file=paste(file.name), sep = "\t" , quote = F, col.names = NA, append = TRUE)
    write(paste("\nPPV:", PPV, "\nNPV:", NPV, "\nTPR:", TPR, "\nTNR:", TNR, "\nMCC:", MCC, "\nACC:", ACC, "\nBACC:", BACC, "\nF1:", F1, "\n\n************************************\n"), file=paste(file.name), append = TRUE)
  }
  
  return(c(BACC, MCC, F1, TPR, TNR))
}


computePValueSDPlot <- function(submission, pValuePredicted, SDPredicted, pValueObserved, complete)
{
  color <- if(complete==TRUE) c(ifelse(pValueObserved<=0.05, "red", "black")) else "blue"
  newData <- data.frame(pValuePredicted, SDPredicted)
  newData <- newData[order(pValuePredicted),]
  simple.name <- gsub(pattern = "_", ".", gsub(pattern = "Submission_", "Submission ", submission))
  # create pch vector
  pch.mut <- as.numeric(fun.eff2)
  pch.mut[as.numeric(fun.eff2) < 1] <- 19
  pch.mut[as.numeric(fun.eff2) > 0] <- 17
  
  if(color[1]!="blue"){
    plot(1:length(newData$pValuePredicted), newData$pValuePredicted, pch=pch.mut, xlab="Predictions", ylab="", ylim=c(0,1), main=simple.name, col=color, cex.axis=0.8, cex.main=0.8)
  } else {
    plot(1:length(newData$pValuePredicted), newData$pValuePredicted, pch=19, xlab="Predictions", ylab="", ylim=c(0,1), main=simple.name, col=color, cex.axis=0.8, cex.main=0.8)
  }
  arrows(1:length(newData$pValuePredicted), newData$pValuePredicted-newData$SDPredicted, 1:length(newData$pValuePredicted), newData$pValuePredicted+newData$SDPredicted, length=0.05, angle=90, code=3, col = "black", lwd=0.3)
}


computePValueSDPlotMultiRow <- function(path.prediction, WT, Experimental.p.value, complete)
{
  pValue.kind <- ifelse(WT == TRUE, "MO_WT", "MO")
  complete.kind <- ifelse(complete == TRUE, "_complete", "")
  png(paste('./results/Figure_S1_ValidationSD_',pValue.kind, complete.kind, '.png', sep=""), width = 7.5, height = 7.5, units = 'in', res = 300)
  par(mar=c(0,2,3,0)+.1, mfrow = c(2, 3), pty = "s", srt = 45) # decide number of plots in each row of the picture and total number of lines
  for(i in 1:length(path.prediction)) {
    CurrentSubmission <- sort(path.prediction)[i]
    submission.name <- gsub("-prediction_file-", ".", gsub("^Group_", "Submission ", gsub(".txt$", "", basename(CurrentSubmission))))
    cat("****************** ", submission.name, " ******************\n")
    
    prediction.table <- read.table(file = CurrentSubmission, sep = "\t", quote = "", header = TRUE)
    prediction.p.value <- if(WT==TRUE) prediction.table$P.value_MO.WT else prediction.table$P.value_MO
    prediction.SD <- if(WT==TRUE) prediction.table$Standard_deviation.1 else prediction.table$Standard_deviation
    
    computePValueSDPlot(submission=submission.name, pValuePredicted=prediction.p.value, SDPredicted=prediction.SD, pValueObserved=Experimental.p.value, complete=complete)
  }
  dev.off()
}

ClassifySubmission <- function(complete.score.matrix)
{
  complete.score.matrix.rank <- data.frame(sapply(1:ncol(complete.score.matrix), function(x) {rank(-complete.score.matrix[, x], ties.method = "average")}))
  rownames(complete.score.matrix.rank) <- c(rownames(complete.score.matrix))
  colnames(complete.score.matrix.rank) <- c(colnames(complete.score.matrix))
  complete.score.matrix.rank <- cbind(complete.score.matrix.rank, rep(NA, nrow(complete.score.matrix.rank)), rep(NA, nrow(complete.score.matrix.rank)))
  colnames(complete.score.matrix.rank) <- c(colnames(complete.score.matrix), "AVG", "Final")
  complete.score.matrix.rank$AVG <- rowMeans(complete.score.matrix.rank, na.rm = TRUE)
  complete.score.matrix.rank$Final <- rank(complete.score.matrix.rank$AVG, ties.method = "average")
  complete.score.matrix.rank.sorted <- complete.score.matrix.rank[order(complete.score.matrix.rank$Final), ]
  write.table(complete.score.matrix.rank.sorted, file='./results/scores_indexes_by_submission_rank.txt', sep = "\t" , quote = F, col.names = NA)
  
  
  complete.score.matrix.sorted <- round(complete.score.matrix[order(complete.score.matrix.rank$Final), ], 2)
  cellcolors<-matrix(NA, ncol(complete.score.matrix.sorted), ncol(complete.score.matrix.sorted))
  threshold = 0
  cellcolors[!is.na(complete.score.matrix.sorted) & complete.score.matrix.sorted > threshold] <- color.scale(complete.score.matrix.sorted[complete.score.matrix.sorted > threshold & !is.na(complete.score.matrix.sorted)], cs1 = 0, cs2 = c(0.9, 0.4), cs3 = 0)
  cellcolors[!is.na(complete.score.matrix.sorted) & complete.score.matrix.sorted < threshold] <- color.scale(complete.score.matrix.sorted[complete.score.matrix.sorted < threshold & !is.na(complete.score.matrix.sorted)], cs1 = c(0.4, 0.9), cs2 = 0, cs3 = 0)
  cellcolors[!is.na(complete.score.matrix.sorted) & complete.score.matrix.sorted == threshold] <- "#FFFFFFFF"
  cellcolors[is.na(complete.score.matrix.sorted)] <- "#000000FF"
  png('./results/Figure_2_scores_indexes_by_submission.png', width = 10, height = 10, units = 'in', res = 300)
  par(mar=c(0,nchar(max(rownames(complete.score.matrix.sorted)))/2+3,nchar(max(colnames(complete.score.matrix.sorted))),0)+.1)
  color2D.matplot(complete.score.matrix.sorted, show.values = 2, axes = FALSE, cellcolors=cellcolors, vcex = 1.6, xlab="", ylab="")
  axis(3,at=0.5:ncol(complete.score.matrix.sorted),labels=colnames(complete.score.matrix.sorted), cex.axis = 1.6)
  axis(2,at=0.5:nrow(complete.score.matrix.sorted),las=2,labels=rev(gsub("_", ".", gsub("Submission_", "Submssion ", rownames(complete.score.matrix.sorted)))), cex.axis = 1.6)
  ppi = 600
  dev.off()
}

ComputeTableRanking <- function(){
  table.score <- read.table("./results/scores_indexes_by_submission.txt", header = T, sep = '\t')
  table.rank <- read.table("./results/scores_indexes_by_submission_rank.txt", header = T, sep ='\t')
  table.rank <- table.rank[order(table.rank$X), ]
  table.sr <- t(sapply(1:nrow(table.rank), function(x){paste(table.rank[x, 2:6 ], " (", table.score[x, 2:6], ")", sep = "")}))
  table.rank[, 2:6] <- table.sr
  table.rank <- table.rank[order(table.rank$Final), ]
  colnames(table.rank)[1] <- "Submission ID"
  write.table(table.rank, file='./results/table_3.txt', sep = "\t" , quote = F, row.names = F)
  
}

ComputeCorrectPredictionsPlot <- function(path.prediction, Experimental.table){
  # create vector with experimental class
  answers.vector <- Experimental.table$Functional.effect_2
  benign <- which(answers.vector=="benign")
  answers.vector[which(answers.vector!="bening")] <- "Disease"
  answers.vector[benign] <- "Benign" 
  # create plot table
  plot.table <- cbind(Experimental.table$Variant, answers.vector)
  for(i in 1:length(path.prediction)) { 
    CurrentSubmission <- sort(path.prediction)[i]
    prediction.table <- read.table(file=CurrentSubmission, sep="\t", header=TRUE, stringsAsFactors = F)
    predicted.class <- prediction.table$Functional_effect
    positive.class <- which(prediction.table$Functional_effect > 0)
    predicted.class <- rep("Benign", length(predicted.class))
    predicted.class[positive.class] <- "Disease"
    plot.table <- cbind(plot.table, predicted.class) 
  }
  plot.table[, 1] <- prediction.table[, 2]
  
  # make table S1
  sub.names <- as.character(c(1.1, 2.1, 4.1, 4.2, 5.1, 6.1))
  matrix.res <-matrix(NA, ncol = 9, nrow = 8)
  matrix.res[1, ] <- c(rep("Submission 4.1", 3), rep("Submission 6.1", 3), rep("Submission 1.1", 3))
  matrix.res[5, ] <- c(rep("Submission 2.1", 3), rep("Submission 4.2", 3), rep("Submission 5.1", 3))
  matrix.res[2, ] <- rep(c(" ", rev(c("Obs. Benign", "Obs. Disease"))), 3)
  matrix.res[6, ] <- rep(c(" ", rev(c("Obs. Benign", "Obs. Disease"))), 3)
  matrix.res[c(2, 3, 4, 6, 7, 8) , c(1, 4, 7)]<- c(" ", rev(c("Pred. Benign", "Pred. Disease")))
  order.vector <- c(3, 6, 1, 2, 4, 5)
  # create confusion matrix
  ConfusionTemp <- function(real.value, predicted.value){
    matr.temp <- matrix(rep(NA, 4), ncol=2)
    matr.temp[1, 1] <- length(intersect(which(real.value=="Disease"), which(predicted.value=="Disease")))
    matr.temp[1, 2] <- length(intersect(which(real.value=="Benign"), which(predicted.value=="Disease")))
    matr.temp[2, 1] <- length(intersect(which(real.value=="Disease"), which(predicted.value=="Benign")))
    matr.temp[2, 2] <- length(intersect(which(real.value=="Benign"), which(predicted.value=="Benign")))
    return(matr.temp)
  }
  matrix.res[c(3, 4), c(8, 9)] <- ConfusionTemp(real.value = answers.vector, predicted.value = plot.table[, 3])
  matrix.res[c(7, 8), c(2, 3)] <- ConfusionTemp(real.value = answers.vector, predicted.value = plot.table[, 4])
  matrix.res[c(3, 4), c(2, 3)] <- ConfusionTemp(real.value = answers.vector, predicted.value = plot.table[, 5])
  matrix.res[c(7, 8), c(5, 6)] <- ConfusionTemp(real.value = answers.vector, predicted.value = plot.table[, 6])
  matrix.res[c(7, 8), c(8, 9)] <- ConfusionTemp(real.value = answers.vector, predicted.value = plot.table[, 7])
  matrix.res[c(3, 4), c(5, 6)] <- ConfusionTemp(real.value = answers.vector, predicted.value = plot.table[, 8])
  write.table(matrix.res, file='./results/Table_S1.txt', sep = "\t" , quote = F, col.names = F, row.names = F)
  
  match.class <- sapply(1:length(answers.vector), function(x){round((length(which(answers.vector[x]==plot.table[x, 3:ncol(plot.table)]))/6)*100, digits = 1)})
  plot.table <- cbind(plot.table[, c(1, 2)], match.class)
  plot.table <- data.frame(plot.table, stringsAsFactors = F)
  colnames(plot.table) <- c("Variant", "Class", "Correct predictions" )
  plot.table$Class <- as.factor(plot.table$Class)
  plot.table$`Correct predictions` <- as.numeric(plot.table$`Correct predictions`)
  plot.table$Variant <- as.factor(plot.table$Variant)
  pdf(file = "./results/Figure_3.pdf", width = 6.5, height = 3.7)
  p <- ggplot(plot.table,aes(x=reorder(Variant, -`Correct predictions`) ,y=`Correct predictions`, group=1)) + 
    geom_line(linetype = "dashed") + geom_point(aes(colour=Class), size=2.5) + theme_gdocs() + 
    theme(axis.text.x = element_text(angle=-90,hjust = 1,vjust = 0.5),legend.position = "top", legend.direction = "horizontal", legend.title = element_blank()) + 
    xlab("") + ylab ("Well predictions [%]") + scale_colour_few()
  print(p)
  dev.off()

  
}

