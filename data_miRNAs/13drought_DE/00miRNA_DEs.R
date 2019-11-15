setwd("/home/zhangt/00projects_2019/01miRNAs_large_scale_maize/data_miRNAs/13_drought_DE")
options(stringsAsFactors = F)

exp_mat <- read.table("../5largematrix/unique_read_RPM.txt", row.names = 1, header = T)
MIRNA_seq <- read.table("MIRNA_isos.txt", header = T)
location <- read.table("miRNAs_and_locations.txt", sep = "\t")
rownames(location) <- location[,2]
final_mirna <- read.table("finalTable.txt", sep = "\t")
final_mirna <- final_mirna[final_mirna[,2]%in%location[,2], ]
rownames(final_mirna) <- final_mirna[,1]
colnames(exp_mat)
sample_information <- read.table("sample_information.txt", sep = "\t")

mirna_num <- final_mirna[,1]
tpm_mat <- t(apply(as.data.frame(mirna_num), 1, function(x){
  seq_out <- MIRNA_seq[MIRNA_seq$pre_name==x, 3]
  exp_out <- exp_mat[seq_out, sample_information[,1]]
  exp_out[is.na(exp_out)] <- 0 
  return(as.numeric(colSums(exp_out)))
}))
rownames(tpm_mat) <- mirna_num
colnames(tpm_mat) <- sample_information[,1]

sample_name <- unique(sample_information[,2])
merge_mat <- matrix(nrow = nrow(tpm_mat), ncol = length(sample_name))
rownames(merge_mat) <- rownames(tpm_mat)
colnames(merge_mat) <- sample_name
for(i in mirna_num){
  value_out <- vector()
  for(j in sample_name){
    value_out <- c(value_out, mean(tpm_mat[i,which(sample_information$V2 %in% j)]))
  }
  merge_mat[i, ] <- value_out
}

# limit_mat_index <- apply(merge_mat, 1, function(x){sum(x>10)>=1})
# merge_mat <- merge_mat[limit_mat_index, ]
# range(apply(merge_mat, 1, max))
# library(pheatmap)
# pheatmap(log2(merge_mat+1), cluster_cols = F)

pvalueCa <- function(x,y){
  if(x==0&y==0){
    res <- 1
  } else if (x==0|y==0){
    res <- 0
  } else {
    res <- 0
    if(x<y){
      for(num in 1:x){
        res=res+log2(y+num)-log2(num)
      }
      res <- 2^(res-x-y-1)
    } else {
      for(num in 1:y){
        res=res+log2(x+num)-log2(num)
      }
      res <- 2^(res-x-y-1)
    }
  }
  res
}

pvalueVec <- vector()
for(i in 1:nrow(merge_mat)){
  resPvalue <- pvalueCa(x=as.numeric(merge_mat[i,1]), y=as.numeric(merge_mat[i,2]))
  pvalueVec <- c(pvalueVec, resPvalue)
}

qvalueVec <- p.adjust(pvalueVec, method ="BH")
logfc <- log2((merge_mat[,2]+1)/(merge_mat[,1]+1))

expValres <- cbind(merge_mat[,1:2], logfc, qvalueVec, "NS")
colnames(expValres)[5] <- "Sig"
B73IndH <- as.numeric(logfc) >= as.numeric(1) & as.numeric(qvalueVec) <= as.numeric(0.05)
B73IndL <- as.numeric(logfc) <= -as.numeric(1) & as.numeric(qvalueVec) <= as.numeric(0.05)
expValres[which(B73IndH),5] <- "Up"
expValres[which(B73IndL),5] <- "Down"
rownames(expValres) <- location[final_mirna[rownames(expValres),2],1]
table(expValres[,5])
write.table(expValres, file = "MIRNAs_DE_result_115.txt", quote = FALSE, sep = "\t", row.names = T)


## based on DEseq
library(corrplot)

count_mat <- read.table("../5largematrix_count/unique_read_count.txt", row.names = 1, header = T)

sample_information <- read.table("sample_information.txt", sep = "\t")
# mirna_count_mat <- matrix(nrow = length(mirna_num), ncol = nrow(sample_information))
# rownames(mirna_count_mat) <- mirna_num
# colnames(mirna_count_mat) <- sample_information[,1]
# for(i in mirna_num){
#   cat(which(mirna_num%in%i), "\n")
#   seq_out <- MIRNA_seq[MIRNA_seq$pre_name==i, 3]
#   exp_out <- count_mat[seq_out, sample_information[,1]]
#   exp_out[is.na(exp_out)] <- 0 
#   mirna_count_mat[i,] <- as.numeric(colSums(exp_out))
# }

count_mat_extract <- t(apply(as.data.frame(mirna_num), 1, function(x){
  seq_out <- MIRNA_seq[MIRNA_seq$pre_name==x, 3]
  exp_out <- count_mat[seq_out, sample_information[,1]]
  exp_out[is.na(exp_out)] <- 0 
  return(as.numeric(colSums(exp_out)))
}))
rownames(count_mat_extract) <- mirna_num
colnames(count_mat_extract) <- sample_information[,1]

count_index <- apply(count_mat_extract, 1, function(x){sum(x>0)>=1})
count_mat_extract <- count_mat_extract[count_index,]

corrplot(cor(count_mat_extract), diag = T)
pheatmap(cor(count_mat_extract))

library(edgeR)
pairs(log2(count_mat_extract+1), pch = ".", lower.panel=NULL)
grp <- factor(sample_information[,2])
# 6483156,3622829,1919373,4859006,9581497,8302080
# d <- DGEList(counts = count_mat_extract, group = grp, lib.size = c(7892605,4516805,4920068,6210636,11479201,10035057))
d <- DGEList(counts = count_mat_extract, group = grp)
d <- calcNormFactors(d, method = "TMM")
head(d$samples)
cols <- as.numeric(d$samples$group)
plotMDS(d, col=cols)

mm <- model.matrix(~-1+grp)
d <- estimateGLMCommonDisp(d,mm)
d <- estimateGLMTrendedDisp(d,mm)
d <- estimateGLMTagwiseDisp(d,mm)
plotBCV(d)
sqrt(d$common.dispersion)
plotMeanVar(d, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)
plotSmear(d, pair=c("Control","Drought"), ylim=c(-5,5))

f <- glmFit(d,mm)
con <- makeContrasts("D-C"=grpDrought-grpControl,levels=colnames(mm))
lrt <- glmLRT(f,contrast=con)
miRNA_DE_result<-topTags(lrt,Inf)$table
miRNA_DE_result<-cbind(miRNAs=rownames(miRNA_DE_result),miRNA_DE_result)

summary(de<-decideTestsDGE(lrt));
detags<-rownames(d)[as.logical(de)];
plotSmear(lrt,de.tags=detags);
abline(h=c(-1,1),col="blue");

# library(DESeq2)
# groupList <- sample_information[,2]
# colData <- data.frame(row.names=colnames(mirna_count_mat), groupList=groupList)
# dds <- DESeqDataSetFromMatrix(countData = mirna_count_mat,colData = colData,design = ~ groupList)
# dds <- DESeq(dds)
# plotDispEsts(dds, main="Dispersion plot")
# res <- results(dds)
# resOrdered <- res[order(res$padj),]
# resOrdered=as.data.frame(resOrdered)

# rld <- rlogTransformation(dds)
# png("RAWvsNORM.png")
# rld <- rlogTransformation(dds)
# exprSet_new=assay(rld)
# par(cex = 0.7)
# n.sample=ncol(mirna_count_mat)
# if(n.sample>40) par(cex = 0.5)
# cols <- rainbow(n.sample*1.2)
# par(mfrow=c(2,2))
# boxplot(mirna_count_mat,  col = cols,main="expression value",las=2)
# boxplot(exprSet_new, col = cols,main="expression value",las=2)
# hist(mirna_count_mat[,1])
# hist(exprSet_new[,1])
# dev.off()
# 
# library(RColorBrewer)
# (mycols <- brewer.pal(8, "Dark2")[1:length(unique(groupList))])
# # Sample distance 
# heatmapsampleDists <- as.matrix(dist(t(exprSet_new)))
# #install.packages("gplots",repos = "http://cran.us.r-project.org")library(gplots)
# png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
# heatmap.2(as.matrix(sampleDists), key=F, 
#           trace="none",col=colorpanel(100, "black", "white"),
#           ColSideColors=mycols[groupList], 
#           RowSideColors=mycols[groupList],margin=c(10, 10), 
#           main="Sample Distance Matrix")
# dev.off()
# png("MA.png")
# DESeq2::plotMA(res, main="DESeq2", ylim=c(-2,2))
# dev.off()
# 
# 
# ### 上面的代码只是为了把6个独立的表达文件给合并成一个表达矩阵
# ## we need to filter the low expression level miRNA
# exprSet=dat[,-1]rownames(exprSet)=dat[,1]
# suppressMessages(library(DESeq2))
# exprSet=ceiling(exprSet)(colData <- data.frame(row.names=colnames(exprSet), groupList=groupList))
# ## DESeq2就是这么简单的用
# dds <- DESeqDataSetFromMatrix(countData = exprSet,colData = colData,design = ~ groupList)
# dds <- DESeq(dds)
# png("qc_dispersions.png", 1000, 1000, pointsize=20)
# plotDispEsts(dds, main="Dispersion plot")
# dev.off()
# res <- results(dds)
# ## 画一些图，相当于做QC吧
# png("RAWvsNORM.png")
# rld <- rlogTransformation(dds)
# exprSet_new=assay(rld)
# par(cex = 0.7)
# n.sample=ncol(exprSet)
# if(n.sample>40) par(cex = 0.5)
# cols <- rainbow(n.sample*1.2)
# par(mfrow=c(2,2))
# boxplot(exprSet,  col = cols,main="expression value",las=2)
# boxplot(exprSet_new, col = cols,main="expression value",las=2)
# hist(exprSet[,1])hist(exprSet_new[,1])
# dev.off()
# library(RColorBrewer)
# (mycols <- brewer.pal(8, "Dark2")[1:length(unique(groupList))])
# # Sample distance 
# heatmapsampleDists <- as.matrix(dist(t(exprSet_new)))
# #install.packages("gplots",repos = "http://cran.us.r-project.org")library(gplots)
# png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
# heatmap.2(as.matrix(sampleDists), key=F, 
#           trace="none",col=colorpanel(100, "black", "white"),
#           ColSideColors=mycols[groupList], 
#           RowSideColors=mycols[groupList],margin=c(10, 10), 
#           main="Sample Distance Matrix")dev.off()
# png("MA.png")
# DESeq2::plotMA(res, main="DESeq2", ylim=c(-2,2))
# dev.off()
# ## 重点就是这里啦，得到了差异分析的结果
# resOrdered <- res[order(res$padj),]
# resOrdered=as.data.frame(resOrdered)
# write.csv(resOrdered,"deseq2.results.csv",quote = F)
# ##下面也是一些图，主要是看看样本之间的差异情况
# library(limma)
# plotMDS(log(counts(dds, normalized=TRUE) + 1))
# plotMDS(log(counts(dds, normalized=TRUE) + 1) - log(t( t(assays(dds)[["mu"]]) / sizeFactors(dds) ) + 1))
# plotMDS( assays(dds)[["counts"]] )  
# ## raw countplotMDS( assays(dds)[["mu"]] ) ##- fitted values.

