# R
setwd("/home/zhangt/00projects_2019/01miRNAs_large_scale_maize/data_miRNAs/12_tissues_RPM/")
options(stringsAsFactors = F)

tukey.biweight <- function(x, c=5, epsilon=0.0001)
{
  m <- median(x)
  s <- median(abs(x - m))
  u <- (x - m) / (c * s + epsilon)
  w <- rep(0, length(x))
  i <- abs(u) <= 1
  w[i] <- ((1 - u^2)^2)[i]
  t.bi <- sum(w * x) / sum(w)
  return(t.bi)
}


exp_ze_mat <- function(final_mirna, samle_information, MIRNA_seq, exp_mat, location){
  out_list <- list()
  mirna_num <- final_mirna[,1]
  tpm_mat <- matrix(nrow = length(mirna_num), ncol = nrow(samle_information))
  rownames(tpm_mat) <- mirna_num
  colnames(tpm_mat) <- samle_information[,1]
  for(i in mirna_num){
    cat(which(mirna_num%in%i), "\n")
    seq_out <- MIRNA_seq[MIRNA_seq$pre_name==i, 3]
    tpm_mat[i,] <- as.numeric(colSums(exp_mat[seq_out, samle_information[,1]]))
  }
  sample_name <- unique(samle_information[,2])
  merge_mat <- matrix(nrow = nrow(tpm_mat), ncol = length(sample_name))
  rownames(merge_mat) <- rownames(tpm_mat)
  colnames(merge_mat) <- sample_name
  for(i in mirna_num){
    value_out <- vector()
    for(j in sample_name){
      value_out <- c(value_out, mean(tpm_mat[i,which(samle_information$V2 %in% j)]))
    }
    merge_mat[i, ] <- value_out
  }
  limit_mat_index <- apply(merge_mat, 1, function(x){sum(x>10)>=1})
  merge_mat <- merge_mat[limit_mat_index, ]
  range(apply(merge_mat, 1, max))
  EXP <- merge_mat
  EXPtukey <- apply(EXP, 1, tukey.biweight)
  
  modEXP <- sapply(1:nrow(EXP), function(aa){return(abs(EXP[aa,] - EXPtukey[aa]))})
  modEXP <- t(modEXP)
  rownames(modEXP) <- rownames(EXP)
  modEXP <- apply(modEXP, 1, unlist)
  modEXP <- t(modEXP)
  
  modH <- apply(modEXP, 1, function(aa){
    bb <- aa/sum(aa)
    sumH <- 0
    for(i in 1:length(bb)){
      sumH <- sumH + bb[i] * log2(bb[i] + 0.00000001)
    }
    return(-sumH)
  })
  EXPZ <- apply(EXP, 1, function(aa){return((aa - mean(aa))/sd(aa))})
  EXPZ <- t(EXPZ)
  EXPZ <- cbind(EXPZ, EXPZ[,1:2])
  EXPZ[,ncol(EXPZ)-1] <- apply(EXPZ[,1:(ncol(EXPZ)-2)], 1, max)
  EXPZ[,ncol(EXPZ)] <- modH
  colnames(EXPZ)[(ncol(EXPZ)-1):ncol(EXPZ)] <- c("maxZ", "modH")
  EXPZ_bc <- EXPZ
  rownames(EXPZ) <- location[final_mirna[rownames(EXPZ),2],1]
  out_list[[1]] <- merge_mat
  out_list[[2]] <- EXPZ
  return(out_list)
}


exp_mat <- read.table("../5largematrix/unique_read_RPM.txt", row.names = 1, header = T)
MIRNA_seq <- read.table("MIRNA_isos.txt", header = T)
location <- read.table("miRNAs_and_locations.txt", sep = "\t")
rownames(location) <- location[,2]
final_mirna <- read.table("finalTable.txt", sep = "\t")
final_mirna <- final_mirna[final_mirna[,2]%in%location[,2], ]
rownames(final_mirna) <- final_mirna[,1]

# merge_tissues
samle_information <- read.table("sample_information_tissues.txt", sep = "\t")
merge_tissues <- exp_ze_mat(final_mirna, samle_information, MIRNA_seq, exp_mat, location)

# split_tissues
# samle_information <- read.table("sample_information.txt", sep = "\t")
# split_tissues <- exp_ze_mat(final_mirna, samle_information, MIRNA_seq, exp_mat, location)

MIRNA <- c("zma-MIR164e", "zma-MIR164h", 
           "zma-MIR169e", "zma-MIR169j", "zma-MIR169o",
           "zma-MIR399g", "zma-MIR399b", "zma-MIR399d",
           "zma-MIR2275a", "zma-MIR2275b", "zma-MIR2275c"
)

EXPZ_tmp1 <- merge_tissues[[2]]
EXPZ_tmp1[intersect(MIRNA, rownames(EXPZ_tmp1)),15:16]
sum(EXPZ_tmp1[,15]>3&EXPZ_tmp1[,16]<1.8)
sort(rownames(EXPZ_tmp1)[EXPZ_tmp1[,15]>3&EXPZ_tmp1[,16]<1.8])
EXPZ_tmp1_result <- EXPZ_tmp1[EXPZ_tmp1[,15]>3&EXPZ_tmp1[,16]<1.8,]
tissues <- apply(EXPZ_tmp1_result[,1:14], 1, function(x){
  colnames(EXPZ_tmp1_result)[which.max(x)]
})
table(tissues)
EXPZ_tmp1_result <- cbind(EXPZ_tmp1_result, tissues)

# EXPZ_tmp2 <- split_tissues[[2]]
# EXPZ_tmp2[intersect(MIRNA, rownames(EXPZ_tmp2)),44:45]
# sum(EXPZ_tmp2[,44]>3.9&EXPZ_tmp2[,45]<3)

EXP <- merge_tissues[[1]]
rownames(EXP) <- location[final_mirna[rownames(EXP),2],1]
EXP <- EXP[EXPZ_tmp1[,15]>3&EXPZ_tmp1[,16]<1.8, ]

write.table(EXP, file = paste0(ncol(EXP), "_samples_",
                               nrow(EXP), "_miRNAs_expssion_values_115.txt"),
            quote = FALSE, sep = "\t", row.names = T)
write.table(EXPZ_tmp1_result, file = paste0(ncol(EXP), "_samples_", 
                                            nrow(EXP), "_miRNAs_shannon_zscore_115.txt"), 
            quote = FALSE, sep = "\t", row.names = T)


pnorm()
qnorm(0.999)
# tukeybiweight <-  function(x, c=5, epsilon=0.0001){
#   list(exprs=apply(x,2,tukey.biweight,c=c,epsilon=epsilon),se.exprs=rep(NA,ncol(x)))
# }
# EXPtukey <- tukeybiweight(t(EXP))


input_data <- read.table("14_samples_94_miRNAs_shannon_zscore_115.txt", sep = "\t", stringsAsFactors = F)
con_select <- sort(unique(input_data[,17]))
con_select <- gsub(" ", ".", con_select)
input_data <- input_data[,con_select]

rownames(input_data) <- gsub("zma-", "", rownames(input_data))
pdf("FigureS3_Heatmap_tissue-specific_miRNA.pdf", height = 10, width = 8)
pheatmap::pheatmap(input_data, cluster_cols = F, border_color = "grey70",
                   color = colorRampPalette(c("white", "#FB1823"))(100),
                   cellwidth = 10, fontsize_row = 7, clustering_method = "ward.D2")
dev.off()
