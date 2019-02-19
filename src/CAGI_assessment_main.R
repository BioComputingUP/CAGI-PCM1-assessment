#!/usr/bin/env/ Rscript

library(ROCR)
library(plotrix)
library(ggthemes)
library(ggplot2)

# set working directory
setwd("../")

# load assessment functions
source("src/validate_format.R")
source("src/Scatterplot_pred_obs.R")
source("src/Calculate_statistics.R")

columnP1 <- 3
columnSD1 <- 4
columnP2 <- 5
columnSD2 <- 6

color.couple <- c('blue', 'yellow')
color.trio <- c("blue", "cornflowerblue", "yellow")

# Import data
exp.val <- read.table('./data/experimental_data/real_proliferation.txt', sep = '\t', header = TRUE, stringsAsFactors = F) #file with experimental data
template.path <- './data/template/template.txt'  # submission template
submission.files <- list.files(path = './data/submissions/', pattern = ".txt", full.names = TRUE)


# prepare experimental table

# remove "<" from experimental values
exp.val[, 3] <- as.numeric(gsub("<", "", exp.val[, 3]))
exp.val[, 4] <- as.numeric(gsub("<", "", exp.val[, 4]))
# create Functional.effect_2 column
fun.eff2 <- exp.val[, 5]
# convert loss of function to 2
fun.eff2[which(fun.eff2 == "loss of function")] <- 2
# convert hypomorph to 1
fun.eff2[which(fun.eff2 == "hypomorph")] <- 1
# convert benign to 0
fun.eff2[which(fun.eff2 == "benign")] <- 0
exp.val <- cbind(exp.val, fun.eff2)
colnames(exp.val)[length(colnames(exp.val))] <- "Functional.effect_2"
exp.val$Functional.effect_2 <- exp.val$Functional.effect
exp.val$Functional.effect <- fun.eff2

#dataset statistics
statistical.data = c()
classes = c("0", "1", "2")
for(class in classes)
  statistical.data = c(statistical.data, length(which(exp.val$Functional.effect==class)))

if(FALSE){
  #STATISTIC PLOTS
  if(FALSE)
  {
    png('./results/database_statistics_threshold_MO_MO_WT.png', width = 10, height = 10, units = 'in', res = 300)
    threshold <- 0.05
    p.value.MO.threshold <- c(length(which(exp.val$p.value.from.MO <= threshold)), length(which(exp.val$p.value.from.MO > threshold)))
    p.value.MO.WT.threshold <- c(length(which(exp.val$p.value.from.MO.WT <= threshold)), length(which(exp.val$p.value.from.MO.WT > threshold)))
    data <- matrix(c(p.value.MO.threshold, p.value.MO.WT.threshold), nrow=2)
    colnames(data) <- c("MO", "MO+WT")
    plt <- barplot(data, col=color.couple, xaxt="n", ylim=c(0, 40), beside=T)
    legend("topright", legend=c(paste("p<=", threshold), paste("p>", threshold)), fill=color.couple, cex=0.8)
    text(mean(plt), (38), labels = "Distribution of p-value", xpd = TRUE, cex=1.6, font=2)
    text(colMeans(plt), -1, labels = colnames(data), xpd = TRUE, cex=1.0)
    dev.off()
  }
  
  png('./results/database_statistics_threshold_MO_.png', width = 3, height = 4, units = 'in', res = 300)
  threshold <- 0.05
  data <- c(length(which(exp.val$p.value.from.MO <= threshold)), length(which(exp.val$p.value.from.MO > threshold)))
  plt <- barplot(data, col=color.couple, xaxt="n", ylim=c(0, 30))
  text(mean(plt), (35), labels = "Distribution of p-value in MO", xpd = TRUE, cex=0.8, font=2)
  text(plt, -2, labels = c("p<=0.05", "p>0.05"), xpd = TRUE, cex=0.8)
  dev.off()
  
  png('./results/database_statistics_threshold_MO_WT_.png', width = 3, height = 4, units = 'in', res = 300)
  threshold <- 0.05
  data <- c(length(which(exp.val$p.value.from.MO.WT <= threshold)), length(which(exp.val$p.value.from.MO.WT > threshold)))
  plt <- barplot(data, col=color.couple, xaxt="n", ylim=c(0, 30))
  text(mean(plt), (35), labels = "Distribution of p-value in MO+WT", xpd = TRUE, cex=0.8, font=2)
  text(plt, -2, labels = c("p<=0.05", "p>0.05"), xpd = TRUE, cex=0.8)
  dev.off()
  
  png('./results/database_statistics_significance_by_MO.png', width = 3, height = 4, units = 'in', res = 300)
  data <- exp.val$p.value.from.MO
  pvalue.001 <- length(which(data<=0.01))
  pvalue.005 <- length(which(data>0.01 & data<=0.05))
  first.column <- ceiling((pvalue.001+pvalue.005)/10)*10
  h <- hist(data, main="P-value distribution for MO", xlab="p-value", ylab="Occurrence", col = color.trio[3], breaks = c(seq(0,1,by=0.05)), freq = TRUE, axes = FALSE, cex.main=0.8, ylim=c(0, first.column), cex.lab=0.6)
  counts=c(pvalue.001+pvalue.005, h$counts[-1])
  axis <- axis(1,at=c(0, 0.05, seq(0.1,1,by=0.1)), cex.axis=0.4, lwd=0.7, labels=c(0, 0.5, seq(0.1,1, by=0.1)))
  axis <- axis(2,at=seq(0,first.column,by=5), cex.axis=0.5, lwd=0.7)
  lw <- h$breaks[1]
  up <- h$breaks[2]
  sz <- h$counts[1]
  rect(lw,0,up,pvalue.001,col=color.trio[1])
  rect(lw,pvalue.001,up,pvalue.001+pvalue.005,col=color.trio[2])
  legend("topright", legend=c("p<=0.01", "0.01<p<=0.05", "p>0.05"), fill=color.trio, cex=0.4)
  dev.off();
  
  png('./results/database_statistics_significance_by_MO_WT.png', width = 3, height = 4, units = 'in', res = 300)
  data <- exp.val$p.value.from.MO.WT
  pvalue.001 <- length(which(data<=0.01))
  pvalue.005 <- length(which(data>0.01 & data<=0.05))
  first.column <- ceiling((pvalue.001+pvalue.005)/10)*10
  h <- hist(data, main="P-value distribution for MO+WT", xlab="p-value", ylab="Occurrence", col = color.trio[3], breaks = c(seq(0,1,by=0.05)), freq = TRUE, axes = FALSE, cex.main=0.8, ylim=c(0, first.column), cex.lab=0.6)
  counts=c(pvalue.001+pvalue.005, h$counts[-1])
  axis <- axis(1,at=c(0, 0.05, seq(0.1,1,by=0.1)), cex.axis=0.4, lwd=0.7, labels=c(0, 0.5, seq(0.1,1, by=0.1)))
  axis <- axis(2,at=seq(0,first.column,by=5), cex.axis=0.5, lwd=0.7)
  lw <- h$breaks[1]
  up <- h$breaks[2]
  sz <- h$counts[1]
  rect(lw,0,up,pvalue.001,col=color.trio[1])
  rect(lw,pvalue.001,up,pvalue.001+pvalue.005,col=color.trio[2])
  legend("topright", legend=c("p<=0.01", "0.01<p<=0.05", "p>0.05"), fill=color.trio, cex=0.4)
  dev.off();
  
    
  jpeg('./results/database_statistics_by_classes.jpg', width = 10, height = 10, units = 'in', res = 300)
  plt <- barplot(statistical.data, col=color.trio, xaxt="n", ylim=c(0, max(statistical.data) + max(statistical.data)/4))
  text(mean(plt), (max(statistical.data) + max(statistical.data)/4), labels = "Distribution of classes", xpd = TRUE, cex=1.6, font=2)
  text(plt, -1, labels = c("Benign", "Hypomorphic", "Loss of Function"), xpd = TRUE, cex=1.4)
  text(plt, statistical.data/2, labels = statistical.data, xpd = TRUE, cex=1.2)
  dev.off()
}

# check submission format
ValidateSubmission(submission.files, nrow=39, ncol=9)

# validate classes in each submission
validateClasses(submission.files)

# Create experimental-observed scatterplots

# P-value MO
ExperimentalVsPredicted(submission.files = submission.files, exp.val = exp.val, column.P = columnP1)
# P-value MO+WT
ExperimentalVsPredicted(submission.files = submission.files, exp.val = exp.val, column.P = columnP2)

#compute basic statistics for predictors
matrices <- ConfusionMatrix(exp.val=exp.val, path.prediction=submission.files)
complete.score.matrix <- as.data.frame(matrices[1])

# P-value MO+WT
CorrelationHeatMap(path.prediction = submission.files)
ClassifySubmission(complete.score.matrix = complete.score.matrix)

# Table 3
ComputeTableRanking()

# Figure 3 
ComputeCorrectPredictionsPlot(path.prediction = submission.files, Experimental.table = exp.val)

#THIS STEP POPULATE AN EMPTY DATAFRAME WITH N RANDOM SUBMISSIONS. EACH ROW IS A MUTATION IN THE ORIGNAL DATASET. THE RATIO OF THE CLASSES IN THE ORIGINAL DATA IS MAINTAINED IN EACH SUBMISSION
if(TRUE){
  number.random.submissions = 10000
  random.submissions <- matrix(NA, nrow = nrow(exp.val), ncol = number.random.submissions)
  for(i in 1:ncol(random.submissions))
    random.submissions[, i] <- sample(c(rep(0, statistical.data[1]), rep(1, statistical.data[2]), rep(2, statistical.data[3])))
  complete.score.matrix.random <- ConfusionMatrixRandomSubmissions(exp.val=exp.val, path.prediction=random.submissions)
  
  png(paste('./results/Figure_S3_Random_BACC_distribution.png', sep=""), width = 6, height = 6, units = 'in', res = 300)
  hist(complete.score.matrix.random$BACC, col="blue", main="Random BACC distribution", xlab = "BACCs", ylab="Occurrence")
  abline(v=0.6647727, col="green", lwd = 2)
  legend("topright", legend = c("Submission 4.1"), fill = c("green"))
  dev.off()
  png(paste('./results/Figure_S3_Random_MCC_distribution.png', sep=""), width = 6, height = 6, units = 'in', res = 300)
  hist(complete.score.matrix.random$MCC, col="blue", main="Random MCC distribution", xlab = "MCCs", ylab="Occurrence")
  abline(v=0.35003330, col="green", lwd = 2)
  legend("topright", legend = c("Submission 4.1"), fill = c("green"))
  dev.off()
  png(paste('./results/Figure_S3_Random_F1_distribution.png', sep=""), width = 6, height = 6, units = 'in', res = 300)
  hist(complete.score.matrix.random$F1, col="blue", main="Random F1 distribution", xlab = "F1s", ylab="Occurrence")
  dev.off()
}

