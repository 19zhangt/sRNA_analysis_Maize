#!/usr/bin/env Rscript
args=commandArgs(T)
library(reshape2)
library(seqinr)

## all expression values
filenames=list.files(path=args[1], full.names=TRUE)
datalist = lapply(filenames, function(x){read.table(file=x,header=F)})
rawTab = do.call(rbind,datalist)
colnames(rawTab) <- c("SRR", "Seq", "RPM")
rawMat <- dcast(rawTab, Seq~SRR, value.var = "RPM")
rownames(rawMat) <- rawMat[,1]
rawMat <- rawMat[,-1]

tmpname <- rownames(rawMat)
if(!is.data.frame(rawMat)){
  cat("\nThis is one sample\n")
  names(rawMat) <- tmpname
}
rawMat[is.na(rawMat)] <- 0
sumval <- apply(rawMat, 1, sum)
leastval <- apply(rawMat, 1, max)
abunIndex <- (sumval>as.numeric(args[3]))|(leastval>as.numeric(args[4]))

rawMat <- rawMat[abunIndex, ]
write.table(rawMat, args[2], quote = F, sep = "\t")

## mireap preparation
mireapMat <- rawMat[nchar(rownames(rawMat))>=as.numeric(args[5])&nchar(rownames(rawMat))<=as.numeric(args[6]), ]
mireapSum <- apply(mireapMat, 1, sum)

write.fasta(sequences = as.list(rownames(mireapMat)), names = paste0("mireap_", seq_len(nrow(mireapMat)), " ", mireapSum), 
            file.out = args[9])

## isomiRs preparation
isoMat <- rawMat[nchar(rownames(rawMat))>=as.numeric(args[7])&nchar(rownames(rawMat))<=as.numeric(args[8]), ]
isoSum <- apply(isoMat, 1, sum)
write.fasta(sequences = as.list(rownames(isoMat)), names = paste0("isoLen_", seq_len(nrow(isoMat)), " ", isoSum), 
            file.out = args[10])
