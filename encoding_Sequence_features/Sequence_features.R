library(BSgenome)
library(Biostrings)
library(dplyr)
library(matrixStats)

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

PSTNP_ForTest <- function(data, n = 3, model = "DC1_Full"){
  kmer <- vector()
  for (i in 1:n){
    kmer <- paste0(rep(kmer, 4), c(rep("A", 4^(i-1)), rep("C", 4^(i-1)), rep("G", 4^(i-1)), rep("U", 4^(i-1))))
  }
  z_score <- readRDS(paste0(model, "_zScoreTable(PSTNP).rds"))
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