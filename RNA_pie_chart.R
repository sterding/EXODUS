library(getopt)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(openxlsx)
library(RColorBrewer)

command <- matrix(
    c("DeRes","D",1,"character",
      "Anno","A",1,"character",
      "Prefix","P",1,"character",
      "outdir","O",2,"character")
      ,byrow=T,ncol=4
)
opt <- getopt(command)

if(!dir.exists(opt$outdir)){
    dir.create(opt$outdir)
}

de <- read.table(opt$DeRes,sep="\t",header=T,stringsAsFactors=F)
anno <- read.xlsx(opt$Anno)

info <- merge(de,anno,by.x="EnsemblID",by.y="ensembl_gene_id",all.x=T)

IGtype <- c("IG_C_gene","IG_D_gene","IG_J_gene","IG_V_gene")
lncRNAtype <- c("3prime_overlapping_ncRNA","antisense_RNA","macro_lncRNA","non_coding","bidirectional_promoter_lncRNA","sense_intronic","sense_overlapping","lincRNA")
ncRNAtype <- c("miRNA","misc_RNA","piRNA","rRNA","siRNA","snRNA","snoRNA","tRNA","vaultRNA","Mt_rRNA","Mt_tRNA")
proteinCodingtype <- c("protein_coding")
Pseudogenetype <- c("IG_pseudogene","IG_C_pseudogene","IG_J_pseudogene","IG_V_pseudogene","polymorphic_pseudogene","processed_pseudogene","pseudogene","transcribed_processed_pseudogene","transcribed_unitary_pseudogene","transcribed_unprocessed_pseudogene","translated_processed_pseudogene","TR_J_pseudogene","TR_V_pseudogene","unitary_pseudogene","unprocessed_pseudogene")
TRtype <- c("TR_C_gene","TR_D_gene","TR_J_gene","TR_V_gene")

typeList <- list(IG_gene=IGtype,lncRNA=lncRNAtype,ncRNA=ncRNAtype,proteinCoding=proteinCodingtype,pseudogene=Pseudogenetype,TR_gene=TRtype)

typeCount <- sapply(typeList,function(x){return(length(which(info$gene_biotype %in% x)))})

others <- nrow(info)-sum(typeCount)

names(others) <- "others"

typeCount <- c(typeCount,others)

typeCount <- typeCount[typeCount!=0]

col <- brewer.pal(length(typeCount), "Set1")

pdf(paste0(opt$outdir,"/",opt$Prefix,"_pieChart.pdf"))

#pie(typeCount,labels=names(typeCount),border="white",col=col)
pie(typeCount,labels=paste(round(typeCount/sum(typeCount),3)*100,"%"),border="white",col=col)

legend("topleft",title="RNA type",names(typeCount),fill=col,cex=0.8,box.lty=0)

dev.off()
