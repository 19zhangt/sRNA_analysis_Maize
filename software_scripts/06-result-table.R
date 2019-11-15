##
library(igraph)
options(stringsAsFactors = F)
# args = commandArgs(trailingOnly=TRUE)

##
mydata <- read.table("data/9mirnas/Conserved3.txt")
g=graph.data.frame(mydata)
newres <- cluster_infomap(g)
names(g[1])
species_name <- names(sort(table(sapply(strsplit(mydata[,1], "-"), function(x)x[1])))[2])

##
genern <- read.table("data/9mirnas/non_Conserved3.txt")
genern[,1] <- sort(genern[,1])
newName <- matrix(nrow = length(names(g[1]))+nrow(genern), ncol = 1)
rownames(newName) <- c(names(g[1]), genern[,1])

rawNum =1
for(i in sort(unique(newres$membership))){
  numindex <- newres$membership%in%i
  geneNam <- names(g[1])[numindex]
  cat(geneNam,"\n")
  if(sum(grepl("MIR", geneNam))>0&!all(grepl("MIR", geneNam))){
    cat(i,"\n")
    newName[geneNam,1] <- geneNam
    novel <- geneNam[!grepl("MIR", geneNam)]
    len_nov <- length(novel)
    gername <- gsub("[a-z]$", "", geneNam[grepl("MIR", geneNam)])[1]
    newName[novel,1] <- unique(paste0(gername, "_N", 1:len_nov))
  } else if(!all(grepl("MIR", geneNam))){
    len_nov <- length(geneNam)
    newName[geneNam,1] <- unique(paste0(species_name, "-MIR_N", rawNum, letters[1:len_nov]))
    rawNum <- rawNum+1
  } else {
    newName[geneNam,1] <- geneNam
  }
}

know_gene <- genern[grepl("MIR", genern[,1]),1]
novel_gene <- genern[grepl("xxx", genern[,1]),1]
newName[know_gene,1] <- know_gene
newName[novel_gene,1] <- paste0(species_name, "-MIR_N", rawNum:(rawNum+length(novel_gene)-1), "a")

write.table(newName, "data/9mirnas/NameChange2.txt", quote = F, sep = "\t", col.names = F)
