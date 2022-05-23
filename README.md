## m1ARegpred: Epitranscriptome Target Prediction of N1-methyladenosine (m1A) Regulators Based on Sequencing Features and Genomic Features

##In this script we demonstrate the method of m1ARegpred framework construction
##Take YTHDC1-Full transcript model as example

##load data
dataset1 <- readRDS("DC1 full transcript dataset 1.rds")
dataset2 <- readRDS("DC1 full transcript dataset 2.rds")
dataset3 <- readRDS("DC1 full transcript dataset 3.rds")
dataset4 <- readRDS("DC1 full transcript dataset 4.rds")

##Four datasets were mixed firstly
library(BSgenome)
require(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
SiteData <- c(dataset1, dataset2, dataset3, dataset4)
SequenceData <- RNAStringSet(DNAStringSet(Views(Hsapiens, SiteData + 20)))
metadata(SequenceData)$label <- c(metadata(dataset1)$label, metadata(dataset2)$label, metadata(dataset3)$label, metadata(dataset4)$label)
SequenceData.index <- list()
SequenceData.index[[1]] <- 1:length(dataset1)
SequenceData.index[[2]] <- (1 + length(dataset1)):(length(dataset1) + length(dataset2))
SequenceData.index[[3]] <- (1 + length(dataset1) + length(dataset2)):(length(dataset1) + length(dataset2) + length(dataset3))
SequenceData.index[[4]] <- (1 + length(dataset1) + length(dataset2) + length(dataset3)):(length(dataset1) + length(dataset2) + length(dataset3) + length(dataset4))

##in this part we showed the sequence encoding methods
library(dplyr)#cbind(), rbind()
library(matrixStats)#about data.frame
library(caret)#a package about machine learning

Kmer <- function(data, k = 2){
  kmer <- vector()
  for (i in 1:k){
    kmer <- paste0(rep(kmer, 4), c(rep("A", 4^(i-1)), rep("C", 4^(i-1)), rep("G", 4^(i-1)), rep("U", 4^(i-1))))
  }
  
  Data <- matrix(ncol = length(kmer), nrow = length(data))
  colnames(Data) <- paste0(kmer, "_freq_Kmer")
  for (i in 1:ncol(Data)){
    Data[, i] <- vcountPattern(kmer[i], data) / (width(data) - k + 1)
  }
  print("Kmer-features finished")
  return(Data)
}

ENAC <- function(data, k = 3){
  Data <- matrix(ncol = 4*(width(data[1]) - k + 1), nrow = length(data))
  colnames(Data) <- paste0("ENAC_pos", rep(1:(ncol(Data)/4), each = 4), "_", rep(c("A", "C", "G", "U"), ncol(Data)/4), "_freq")
  for (i in 1:(ncol(Data)/4)){
    Data[, 4*i - 3] <- vcountPattern("A", subseq(data, start = i, width = k)) / k
    Data[, 4*i - 2] <- vcountPattern("C", subseq(data, start = i, width = k)) / k
    Data[, 4*i - 1] <- vcountPattern("G", subseq(data, start = i, width = k)) / k
    Data[, 4*i] <- vcountPattern("U", subseq(data, start = i, width = k)) / k
  }
  print("ENAC-features finished")
  return(Data)
}

PS <- function(data, n = 2){
  mc <- vector()
  for (i in 1:n){
    mc <- paste0(rep(mc, 4), c(rep("A", 4^(i-1)), rep("C", 4^(i-1)), rep("G", 4^(i-1)), rep("U", 4^(i-1))))
  }
  
  Data <- matrix(ncol = length(mc) * (width(data[1]) - n + 1), nrow = length(data))
  colnames(Data) <- paste0("pos", rep(1:(width(data[1]) - n + 1), each = length(mc)), " = ", rep(mc, (width(data[1]) - n + 1)))
  for (i in 1:ncol(Data)){
    Data[, i] <- as.numeric(as.character(subseq(data, start = ((i - 1)%/%length(mc) + 1 ), width = n)) == rep(mc, (width(data[1]) - n + 1))[i])
  }
  print("PS-features finished")
  return(Data)
}

PseE2P <- function(data, k = 3){
  kmer <- vector()
  for (i in 1:k){
    kmer <- paste0(rep(kmer, 4), c(rep("A", 4^(i-1)), rep("C", 4^(i-1)), rep("G", 4^(i-1)), rep("U", 4^(i-1))))
  }
  freq <- matrix(ncol = length(kmer), nrow = length(data))
  for (i in 1:ncol(freq)){
    freq[, i] <- vcountPattern(kmer[i], data) / (width(data) - k + 1)
  }
  
  E2P_list <- rep(0, length(kmer))
  names(E2P_list) <- kmer
  for (i in 1:length(E2P_list)){
    E2P_list[i] <- 0.1260 * sum(unlist(as.vector(strsplit(names(E2P_list)[i],"", fixed = TRUE))) == "A") + 0.1340 * sum(unlist(as.vector(strsplit(names(E2P_list)[i],"", fixed = TRUE))) == "C") + 0.0806 * sum(unlist(as.vector(strsplit(names(E2P_list)[i],"", fixed = TRUE))) == "G") + 0.1335 * sum(unlist(as.vector(strsplit(names(E2P_list)[i],"", fixed = TRUE))) == "U")
  }
  E2Pvalue <-t(replicate(length(data), E2P_list))
  
  Data <- freq * E2Pvalue
  colnames(Data) <- paste0(names(E2P_list), "_E2P*freq")
  print("PseE2P-features finished")
  return(Data)
}

DPCP <- function(data){
  kmer <- vector()
  for (i in 1:2){
    kmer <- paste0(rep(kmer, 4), c(rep("A", 4^(i-1)), rep("C", 4^(i-1)), rep("G", 4^(i-1)), rep("U", 4^(i-1))))
  }
  freq <- matrix(ncol = 9 * length(kmer), nrow = length(data))
  for (i in 1:length(kmer)){
    freq[, c((9 * i - 8):(9 * i))] <- vcountPattern(kmer[i], data) / (width(data) - 1)
  }
  
  PCvalue_table <- read.csv("PC value.csv")
  PCvalue <- matrix(nrow = length(data), ncol = (ncol(PCvalue_table) - 1) * nrow(PCvalue_table))
  for (i in 1:(ncol(PCvalue_table) - 1)){
    PCvalue[, c((9 * i - 8):(9 * i))] <- matrix(rep(PCvalue_table[, (i + 1)], each = length(data)), nrow = length(data))
  }
  
  Data <- freq * PCvalue
  colnames(Data) <- paste0("DPCP_", rep(kmer, each = 9), "_", rep(PCvalue_table[, 1], length(kmer)), "*freq")
  print("DPCP-features finished")
  return(Data)
}

NM1 <- function(data){
  data <- as.character(data)
  NucleoMatrixfirst <- unlist(as.vector(strsplit(data[1],"",fixed = TRUE)))
  NucleoMatrix <- matrix(NA,nrow = length(data),ncol = length(NucleoMatrixfirst))
  for(i in 1:length(data) ){
    NucleoMatrix[i,] <- unlist(as.vector(strsplit(data[i],"",fixed = TRUE)))
  }
  return(NucleoMatrix)
}

NCP <- function(data){
  MA <- NM1(data)
  MA2 <- matrix(NA,ncol =(3*ncol(MA)),nrow = nrow(MA))
  colnames(MA2) <- paste0("cp", 1:ncol(MA2))
  for (i in 1:ncol(MA)){
    for(j in 1:nrow(MA)){
      if (MA[j, i] == "A"){
        MA2[j, (i*3-2):(i*3)] = c(1, 1, 1)
      }else if (MA[j, i] == "C"){
        MA2[j, (i*3-2):(i*3)] = c(0, 1, 0)
      }else if (MA[j, i] == "G"){
        MA2[j, (i*3-2):(i*3)] = c(1, 0, 0)
      }else if (MA[j, i] == "T" | MA[j, i] == "U"){
        MA2[j, (i*3-2):(i*3)] = c(0, 0, 1)
      }else{
        MA2[j, (i*3-2):(i*3)] = c(0, 0, 0)
      }
    }
  }
  print("NCP-features finished")
  return(MA2)
}

ANP <- function(Data,NTYPE="RNA") {
  sequences_M <- NM1(Data)
  if(NTYPE=="RNA"){
    U="U"
  }else{
    U="T"
  }
  N = ncol(sequences_M)
  
  cumFreq_A <- rowCumsums(matrix(as.numeric( sequences_M == "A" ), ncol = N, byrow = F))
  cumFreq_T <- rowCumsums(matrix(as.numeric( sequences_M ==  U ), ncol = N, byrow = F))
  cumFreq_C <- rowCumsums(matrix(as.numeric( sequences_M == "C" ), ncol = N, byrow = F))
  cumFreq_G <- rowCumsums(matrix(as.numeric( sequences_M == "G" ), ncol = N, byrow = F))
  
  cumFreq_combined <- matrix(0,ncol = N, nrow = length(Data))
  cumFreq_combined[sequences_M == "A"] <- cumFreq_A[sequences_M == "A"]
  cumFreq_combined[sequences_M == U] <- cumFreq_T[sequences_M == U]
  cumFreq_combined[sequences_M == "C"] <- cumFreq_C[sequences_M == "C"]
  cumFreq_combined[sequences_M == "G"] <- cumFreq_G[sequences_M == "G"]
  
  cumFreq_combined <- t(t(cumFreq_combined) / seq_len(N))
  colnames(cumFreq_combined ) <-  paste0("cumFreq_",seq_len(N))
  print("ANP-features finished")
  return(cumFreq_combined)
}

PSTNP <- function(data, n = 3){
  kmer <- vector()
  for (i in 1:n){
    kmer <- paste0(rep(kmer, 4), c(rep("A", 4^(i-1)), rep("C", 4^(i-1)), rep("G", 4^(i-1)), rep("U", 4^(i-1))))
  }
  z_score <- matrix(nrow = width(data[1]) - n + 1, ncol = length(kmer))
  colnames(z_score) <- kmer
  for (i in 1:nrow(z_score)){
    x <- as.character(subseq(data, start = i, width = n))
    for (j in 1:ncol(z_score))
      z_score[i, j] <- sum(x[which(metadata(data)$label == "+")] == kmer[j]) / length(which(metadata(data)$label == "+")) - sum(x[which(metadata(data)$label == "-")] == kmer[j]) / length(which(metadata(data)$label == "-"))
  }
  
  Data <- matrix(ncol = width(data[1]) - n + 1, nrow = length(data))
  colnames(Data) <- paste0("PSTNP_pos", c(1:(width(data[1]) - n + 1)))
  for (i in 1:ncol(Data)){
    Data[, i] <- z_score[i, ][as.character(subseq(data, start = i, width = n))]
  }
  print("PSTNP-features finished")
  return(Data)
}

PseNAC <- function(data, n = 5){
  kmer <- vector()
  for (i in 1:2){
    kmer <- paste0(rep(kmer, 4), c(rep("A", 4^(i-1)), rep("C", 4^(i-1)), rep("G", 4^(i-1)), rep("U", 4^(i-1))))
  }
  PCvalue_table <- read.csv("PC value.csv")
  NormalizedPCvalue_table <- PCvalue_table
  for (i in 1:nrow(NormalizedPCvalue_table)){
    NormalizedPCvalue_table[i, 2:17] <- (PCvalue_table[i, 2:17] - mean(as.numeric(PCvalue_table[i, 2:17]))) / sd(as.numeric(PCvalue_table[i, 2:17]))
  }
  
  CorrelationTable <- vector(length = 16 * 16)
  name <- rep(kmer, each = 16)
  for (i in 1:16){
    names(CorrelationTable)[(16*i - 15) : (16*i)] <- paste0(name[(16*i - 15) : (16*i)], "-", kmer)
  }
  for (i in 1:16){
    for (j in 1:16){
      CorrelationTable[16 * i - 16 + j] <- mean((NormalizedPCvalue_table[1:6, kmer[i]] - NormalizedPCvalue_table[1:6, kmer[j]]) ^ 2)
    }
  }
  
  Data <- matrix(ncol = n, nrow = length(data))
  colnames(Data) <- paste0(1:n, "-tier correlation")
  theta <- matrix(nrow = length(data), ncol = (width(data)[[1]] - n - 1))
  for (i in 1:n){
    for (j in 1:(width(data)[[1]] - n - 1)){
      theta[, j] <- CorrelationTable[paste0(as.character(subseq(data, start = j, width = 2)), "-", as.character(subseq(data, start = j + n, width = 2)))]
    }
    Data[, i] <- rowMeans(theta)
  }
  print("PseNAC-features finished")
  return(Data)
}

Sequence_Features <- cbind(Kmer(SequenceData), ENAC(SequenceData), PS(SequenceData), PseE2P(SequenceData), DPCP(SequenceData), NCP(SequenceData), ANP(SequenceData), PSTNP(SequenceData), PseNAC(SequenceData))


##In this part we show the method of generating genome features
##it is done by the function of the m6ALogisticModel package
##these packages are required for the function
library(m6ALogisticModel)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)

matureSE <- SummarizedExperiment()
rowRanges(matureSE) <- SiteData
Additional_features_hg19 = list(
  HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
  YTHDC1_TREW = YTHDC1_TREW_gr,
  YTHDF1_TREW = YTHDF1_TREW_gr,
  YTHDF2_TREW = YTHDF2_TREW_gr,
  miR_targeted_genes = miR_targeted_genes_grl,
  TargetScan = TargetScan_hg19_gr,
  Verified_miRtargets = verified_targets_gr,
  METTL3_TREW = METTL3_TREW,
  METTL14_TREW = METTL14_TREW,
  WTAP_TREW = WTAP_TREW,
  METTL16_CLIP = METTL16_CLIP,
  ALKBH5_PARCLIP = ALKBH5_PARCLIP,
  FTO_CLIP = FTO_CLIP,
  FTO_eCLIP = FTO_eCLIP
)
matureFE <- predictors_annot(se = matureSE,
                             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                             bsgnm = Hsapiens,
                             fc = fitCons.UCSC.hg19,
                             pc = phastCons100way.UCSC.hg19,
                             struct_hybridize = Struc_hg19,
                             feature_lst = Additional_features_hg19,
                             hk_genes_list = HK_hg19_eids,
                             genes_ambiguity_method = "average")
Genomic_Features_c <- data.matrix(mcols(matureFE))
Genomic_Features <- Genomic_Features_c[, c(1:41, 46:51, 56, 66:73)]##select the needed genome features


##In this part we show the method of generating training data, validation data, and testing data
Features_c <- cbind(Sequence_Features, Genomic_Features)
Features <- Features_c[!duplicated(SiteData), ]#remove the repeat site in mix-sample data
Features_labels <- metadata(SequenceData)$label[!duplicated(SiteData)]
set.seed(204)
inTrain <- createDataPartition(Features_labels, p = .8, list = FALSE)##1/5 of the data were selected as the testing data
trainx <- Features[inTrain, ] %>% as.data.frame
testx <- Features[-inTrain, ] %>% as.data.frame

trainy_c = Features_labels[inTrain]
trainy <- vector(length = length(trainy_c))
trainy[which((trainy_c == "+"))] <- "Bind"; trainy[which((trainy_c == "-"))] <- "NoBind"

testy_c = Features_labels[-inTrain]
testy <- vector(length = length(testy_c))
testy[which((testy_c == "+"))] <- 1; testy[which((testy_c == "-"))] <- 0


##Then 1/5 of the remaining data were selected as the validation data
set.seed(666)
inValidation = createDataPartition(trainy, p = 0.8, list = FALSE)
trainset <- trainx[inValidation, ]
trainlabel <- trainy[inValidation]

validationset <- trainx[-inValidation, ]
validationlabel <- vector(length = length(trainy[-inValidation]))
validationlabel[which(trainy[-inValidation] == "Bind")] <- 1; validationlabel[which(trainy[-inValidation] == "NoBind")] <- 0


##In this part we showed the function of Fscore
F_score <- function(data, label){
  F_score <- vector(length = ncol(data))
  names(F_score) <- colnames(data)
  for (i in 1:ncol(data)){
    x_positive <- data[label == "Bind", i]
    x_negative <- data[label == "NoBind", i]
    x <- data[, i]
    F_score[i] <- ((mean(x_positive) - mean(x))^2 +(mean(x_negative) - mean(x))^2) / (var(x_positive) + var(x_negative))
  }
  return(F_score)
}
F_importance <- sort(F_score(trainset[, -(877:1020)], trainlabel), decreasing = T)
selectedtrainset <- trainset[, names(F_importance[1:30])]#select the Top 30 features
selectedtrainset$label <- factor(trainlabel, labels=c("Bind","NoBind"))


##In this part we showed the method of feature selection
library(dplyr)
library(caret)#a package of machine learning
library(ROCR)#for calculating AUC score and AUPRC score
auc <- vector(length = 30)
aucpr <- vector(length = 30)
names(auc) <- c(1:30); names(aucpr) <- c(1:30)
for (i in 1:30){
  trainingset <- selectedtrainset[, c(1:i, 31)]
  set.seed(200)
  SVM <- train(label ~ ., data = trainingset, 
               method = "svmRadial", 
               preProc = c("center", "scale"),
               trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
  predictedProbs <- predict(SVM, newdata = validationset, type = "prob")
  pred <- prediction(predictedProbs$Bind, validationlabel)
  auc[i] <- attr(performance(pred,"auc"), "y.values")[[1]]
  aucpr[i] <- attr(performance(pred,"aucpr"), "y.values")[[1]]
  print(i)
}


##In this part we show the process of evaluating the model using independent test set
feature_number <- 26##the feature number were selected when the score stop growing
trainingset <- selectedtrainset[, 1:feature_number]
trainingset$label <- selectedtrainset$label
set.seed(7)
SVM <- train(label ~ ., data = trainingset, 
             method = "svmRadial", 
             preProc = c("center", "scale"),
             trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))

predictedProbs <- predict(SVM, newdata = testx, type = "prob")
pred <- prediction(predictedProbs$Bind, testy)
CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
CM <- confusionMatrix(predict(SVM, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
result <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])
perf1 = performance(pred, measure="tpr", x.measure="fpr" ) 
plot(perf1)
perf2 <- performance(pred, "prec", "rec")
plot(perf2)


##build the model using conventional encoding methods for comparison
e <- list()
e[[1]] <- trainset[, 1:16]
e[[2]] <- trainset[, 813:876]
e[[3]] <- trainset[, 1144:1184]
e[[4]] <- trainset[, 1185:1223]
convention <- matrix(nrow = 4, ncol = 5)
rownames(convention) <- c("Kmer", "E2P", "ANP", "PSTNP")
colnames(convention) <- c("AUC", "AUPRC", "Accuracy", "Sensitivity", "Specificity")
for (i in 1:4){
  a <- e[[i]]
  nzv1 <- nearZeroVar(a, saveMetrics= TRUE)
  x <- a[, !nzv1$nzv]
  x$label <- selectedtrainset$label
  set.seed(1200)
  Comparision <- train(label ~ ., data = x, 
                       method = "svmRadial", 
                       preProc = c("center", "scale"),
                       trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
  predictedProbs <- predict(Comparision, newdata = testx, type = "prob")
  pred <- prediction(predictedProbs$Bind, testy)
  CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
  CM <- confusionMatrix(predict(Comparision, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
  convention[i, ] <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])
  print(i)
}


##compare different machine learning algorithms
trainingset <- selectedtrainset[, 1:feature_number]
trainingset$label <- selectedtrainset$label
model <- matrix(nrow = 4, ncol = 5)
rownames(model) <- c("gbm", "baye", "RF", "kknn")
colnames(model) <- c("AUC", "AUPRC", "Accuracy", "Sensitivity", "Specificity")

set.seed(727)
gbm <- train(label ~ ., data = trainingset, 
             method = "gbm", 
             preProc = c("center", "scale"),
             trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
predictedProbs <- predict(gbm, newdata = testx, type = "prob")
pred <- prediction(predictedProbs$Bind, testy)
CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
CM <- confusionMatrix(predict(gbm, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
model[1, ] <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])

set.seed(727)
bayes <- train(label ~ ., data = trainingset, 
               method = "naive_bayes", 
               preProc = c("center", "scale"),
               trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
predictedProbs <- predict(bayes, newdata = testx, type = "prob")
pred <- prediction(predictedProbs$Bind, testy)
CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
CM <- confusionMatrix(predict(bayes, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
model[2, ] <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])

set.seed(727)
RF <- train(label ~ ., data = trainingset, 
            method = "rf", 
            preProc = c("center", "scale"),
            trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
predictedProbs <- predict(RF, newdata = testx, type = "prob")
pred <- prediction(predictedProbs$Bind, testy)
CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
CM <- confusionMatrix(predict(RF, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
model[3, ] <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])

set.seed(727)
kknn <- train(label ~ ., data = trainingset, 
              method = "kknn", 
              preProc = c("center", "scale"),
              trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
predictedProbs <- predict(kknn, newdata = testx, type = "prob")
pred <- prediction(predictedProbs$Bind, testy)
CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
CM <- confusionMatrix(predict(kknn, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
model[4, ] <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])


##In this part we show the evaluation process of cross-sample testing
##The first data set as testing set
cross_sample <- list()
cross_sample[[1]] <- matrix(nrow = 5, ncol = 5)
colnames(cross_sample[[1]]) <- c("AUC", "AUPRC", "Accuracy", "Sensitivity", "Specificity")
rownames(cross_sample[[1]]) <- c("dataset1-m1a", "dataset1-Kmer", "dataset1-E2P", "dataset1-ANP", "dataset1-PSTNP")

Features_c <- cbind(Sequence_Features, Genomic_Features)[c(SequenceData.index[[2]], SequenceData.index[[3]], SequenceData.index[[4]]), ]
Features <- Features_c[!duplicated(SiteData[c(SequenceData.index[[2]], SequenceData.index[[3]], SequenceData.index[[4]])]), ]
Features_labels_c <- metadata(SequenceData)$label[c(SequenceData.index[[2]], SequenceData.index[[3]], SequenceData.index[[4]])][!duplicated(SiteData[c(SequenceData.index[[2]], SequenceData.index[[3]], SequenceData.index[[4]])])]
Features_labels <- vector(length = length(Features_labels_c))
Features_labels[which((Features_labels_c == "+"))] <- "Bind"; Features_labels[which((Features_labels_c == "-"))] <- "NoBind"
testx <- cbind(Sequence_Features, Genomic_Features)[SequenceData.index[[1]], ]
testy <- c(rep(1, sum(metadata(SequenceData)$label[SequenceData.index[[1]]] == "+")), rep(0, sum(metadata(SequenceData)$label[SequenceData.index[[1]]] == "-")))

trainingset <- Features[, names(F_importance[1:feature_number])] %>% as.data.frame
trainingset$label <- Features_labels
set.seed(7)
SVM <- train(label ~ ., data = trainingset, 
             method = "svmRadial", 
             preProc = c("center", "scale"),
             trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
predictedProbs <- predict(SVM, newdata = testx, type = "prob")
pred <- prediction(predictedProbs$Bind, testy)
CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
CM <- confusionMatrix(predict(SVM, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
result <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])
cross_sample[[1]][1, ] <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])

e <- list()
e[[1]] <- Features[, 1:16]
e[[2]] <- Features[, 813:876]
e[[3]] <- Features[, 1144:1184]
e[[4]] <- Features[, 1185:1223]
for (i in 1:4){
  a <- e[[i]] %>% as.data.frame
  nzv1 <- nearZeroVar(a, saveMetrics= TRUE)
  x <- a[, !nzv1$nzv]
  x$label <- Features_labels
  set.seed(226)
  Comparision <- train(label ~ ., data = x, 
                       method = "svmRadial", 
                       preProc = c("center", "scale"),
                       trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
  predictedProbs <- predict(Comparision, newdata = testx, type = "prob")
  pred <- prediction(predictedProbs$Bind, testy)
  CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
  CM <- confusionMatrix(predict(Comparision, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
  cross_sample[[1]][(i + 1), ]<- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])
  print(i)
}


##The second data set as testing set
cross_sample[[2]] <- matrix(nrow = 5, ncol = 5)
colnames(cross_sample[[2]]) <- c("AUC", "AUPRC", "Accuracy", "Sensitivity", "Specificity")
rownames(cross_sample[[2]]) <- c("dataset2-m1a", "dataset2-Kmer", "dataset2-E2P", "dataset2-ANP", "dataset2-PSTNP")

Features_c <- cbind(Sequence_Features, Genomic_Features)[c(SequenceData.index[[1]], SequenceData.index[[3]], SequenceData.index[[4]]), ]
Features <- Features_c[!duplicated(SiteData[c(SequenceData.index[[1]], SequenceData.index[[3]], SequenceData.index[[4]])]), ]
Features_labels_c <- metadata(SequenceData)$label[c(SequenceData.index[[1]], SequenceData.index[[3]], SequenceData.index[[4]])][!duplicated(SiteData[c(SequenceData.index[[1]], SequenceData.index[[3]], SequenceData.index[[4]])])]
Features_labels <- vector(length = length(Features_labels_c))
Features_labels[which((Features_labels_c == "+"))] <- "Bind"; Features_labels[which((Features_labels_c == "-"))] <- "NoBind"
testx <- cbind(Sequence_Features, Genomic_Features)[SequenceData.index[[2]], ]
testy <- c(rep(1, sum(metadata(SequenceData)$label[SequenceData.index[[2]]] == "+")), rep(0, sum(metadata(SequenceData)$label[SequenceData.index[[2]]] == "-")))

trainingset <- Features[, names(F_importance[1:feature_number])] %>% as.data.frame
trainingset$label <- Features_labels
set.seed(7)
SVM <- train(label ~ ., data = trainingset, 
             method = "svmRadial", 
             preProc = c("center", "scale"),
             trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
predictedProbs <- predict(SVM, newdata = testx, type = "prob")
pred <- prediction(predictedProbs$Bind, testy)
CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
CM <- confusionMatrix(predict(SVM, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
result <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])
cross_sample[[2]][1, ] <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])

e <- list()
e[[1]] <- Features[, 1:16]
e[[2]] <- Features[, 813:876]
e[[3]] <- Features[, 1144:1184]
e[[4]] <- Features[, 1185:1223]
for (i in 1:4){
  a <- e[[i]] %>% as.data.frame
  nzv1 <- nearZeroVar(a, saveMetrics= TRUE)
  x <- a[, !nzv1$nzv]
  x$label <- Features_labels
  set.seed(226)
  Comparision <- train(label ~ ., data = x, 
                       method = "svmRadial", 
                       preProc = c("center", "scale"),
                       trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
  predictedProbs <- predict(Comparision, newdata = testx, type = "prob")
  pred <- prediction(predictedProbs$Bind, testy)
  CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
  CM <- confusionMatrix(predict(Comparision, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
  cross_sample[[2]][(i + 1), ]<- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])
  print(i)
}


##The third data set as testing set
cross_sample[[3]] <- matrix(nrow = 5, ncol = 5)
colnames(cross_sample[[3]]) <- c("AUC", "AUPRC", "Accuracy", "Sensitivity", "Specificity")
rownames(cross_sample[[3]]) <- c("dataset3-m1a", "dataset3-Kmer", "dataset3-E2P", "dataset3-ANP", "dataset3-PSTNP")

Features_c <- cbind(Sequence_Features, Genomic_Features)[c(SequenceData.index[[1]], SequenceData.index[[2]], SequenceData.index[[4]]), ]
Features <- Features_c[!duplicated(SiteData[c(SequenceData.index[[1]], SequenceData.index[[2]], SequenceData.index[[4]])]), ]
Features_labels_c <- metadata(SequenceData)$label[c(SequenceData.index[[1]], SequenceData.index[[2]], SequenceData.index[[4]])][!duplicated(SiteData[c(SequenceData.index[[1]], SequenceData.index[[2]], SequenceData.index[[4]])])]
Features_labels <- vector(length = length(Features_labels_c))
Features_labels[which((Features_labels_c == "+"))] <- "Bind"; Features_labels[which((Features_labels_c == "-"))] <- "NoBind"
testx <- cbind(Sequence_Features, Genomic_Features)[SequenceData.index[[3]], ]
testy <- c(rep(1, sum(metadata(SequenceData)$label[SequenceData.index[[3]]] == "+")), rep(0, sum(metadata(SequenceData)$label[SequenceData.index[[3]]] == "-")))

trainingset <- Features[, names(F_importance[1:feature_number])] %>% as.data.frame
trainingset$label <- Features_labels
set.seed(7)
SVM <- train(label ~ ., data = trainingset, 
             method = "svmRadial", 
             preProc = c("center", "scale"),
             trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
predictedProbs <- predict(SVM, newdata = testx, type = "prob")
pred <- prediction(predictedProbs$Bind, testy)
CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
CM <- confusionMatrix(predict(SVM, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
result <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])
cross_sample[[3]][1, ] <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])

e <- list()
e[[1]] <- Features[, 1:16]
e[[2]] <- Features[, 813:876]
e[[3]] <- Features[, 1144:1184]
e[[4]] <- Features[, 1185:1223]
for (i in 1:4){
  a <- e[[i]] %>% as.data.frame
  nzv1 <- nearZeroVar(a, saveMetrics= TRUE)
  x <- a[, !nzv1$nzv]
  x$label <- Features_labels
  set.seed(226)
  Comparision <- train(label ~ ., data = x, 
                       method = "svmRadial", 
                       preProc = c("center", "scale"),
                       trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
  predictedProbs <- predict(Comparision, newdata = testx, type = "prob")
  pred <- prediction(predictedProbs$Bind, testy)
  CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
  CM <- confusionMatrix(predict(Comparision, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
  cross_sample[[3]][(i + 1), ]<- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])
  print(i)
}


##The fourth data set as testing set
cross_sample[[4]] <- matrix(nrow = 5, ncol = 5)
colnames(cross_sample[[4]]) <- c("AUC", "AUPRC", "Accuracy", "Sensitivity", "Specificity")
rownames(cross_sample[[4]]) <- c("dataset4-m1a", "dataset4-Kmer", "dataset4-E2P", "dataset4-ANP", "dataset4-PSTNP")

Features_c <- cbind(Sequence_Features, Genomic_Features)[c(SequenceData.index[[1]], SequenceData.index[[2]], SequenceData.index[[3]]), ]
Features <- Features_c[!duplicated(SiteData[c(SequenceData.index[[1]], SequenceData.index[[2]], SequenceData.index[[3]])]), ]
Features_labels_c <- metadata(SequenceData)$label[c(SequenceData.index[[1]], SequenceData.index[[2]], SequenceData.index[[3]])][!duplicated(SiteData[c(SequenceData.index[[1]], SequenceData.index[[2]], SequenceData.index[[3]])])]
Features_labels <- vector(length = length(Features_labels_c))
Features_labels[which((Features_labels_c == "+"))] <- "Bind"; Features_labels[which((Features_labels_c == "-"))] <- "NoBind"
testx <- cbind(Sequence_Features, Genomic_Features)[SequenceData.index[[4]], ]
testy <- c(rep(1, sum(metadata(SequenceData)$label[SequenceData.index[[4]]] == "+")), rep(0, sum(metadata(SequenceData)$label[SequenceData.index[[4]]] == "-")))

trainingset <- Features[, names(F_importance[1:feature_number])] %>% as.data.frame
trainingset$label <- Features_labels
set.seed(7)
SVM <- train(label ~ ., data = trainingset, 
             method = "svmRadial", 
             preProc = c("center", "scale"),
             trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
predictedProbs <- predict(SVM, newdata = testx, type = "prob")
pred <- prediction(predictedProbs$Bind, testy)
CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
CM <- confusionMatrix(predict(SVM, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
result <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])
cross_sample[[4]][1, ] <- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])

e <- list()
e[[1]] <- Features[, 1:16]
e[[2]] <- Features[, 813:876]
e[[3]] <- Features[, 1144:1184]
e[[4]] <- Features[, 1185:1223]
for (i in 1:4){
  a <- e[[i]] %>% as.data.frame
  nzv1 <- nearZeroVar(a, saveMetrics= TRUE)
  x <- a[, !nzv1$nzv]
  x$label <- Features_labels
  set.seed(226)
  Comparision <- train(label ~ ., data = x, 
                       method = "svmRadial", 
                       preProc = c("center", "scale"),
                       trControl = trainControl(method = "cv", number = 5, classProbs =  TRUE))
  predictedProbs <- predict(Comparision, newdata = testx, type = "prob")
  pred <- prediction(predictedProbs$Bind, testy)
  CM_label <- vector(length = length(testy)); CM_label[which(testy == 1)] <- "Bind"; CM_label[which(testy == 0)] <- "NoBind"
  CM <- confusionMatrix(predict(Comparision, newdata = testx), factor(CM_label, labels=c("Bind","NoBind")))
  cross_sample[[4]][(i + 1), ]<- c(attr(performance(pred,"auc"), "y.values")[[1]], attr(performance(pred,"aucpr"), "y.values")[[1]], CM$overall[1], CM$byClass[1:2])
  print(i)
}
