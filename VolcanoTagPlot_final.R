library(EnhancedVolcano)
library(airway)
library(magrittr)
library(Cairo, quietly=TRUE)
file <- commandArgs(T)
res1 <- read.table(file[1],sep="\t",header=T,row.names=1,stringsAsFactors=F)
markers <- read.table(file[2],sep="\t",header=T,stringsAsFactors=F)
res1$anno <- NA
res1[markers$ENSEMBL.ID,]$anno <- markers$GENE.SYMBOL
p <- EnhancedVolcano(res1,lab=res1$anno,selectLab=markers$GENE.SYMBOL,x="log2FoldChange",y="pvalue",xlim=c(-10,10),ylim=c(0,5),pCutoff=0.05,FCcutoff=1,colAlpha=1,DrawConnectors = TRUE,widthConnectors=0.1,colConnectors='grey50')
ggsave(p,file=file[3])
