library(getopt)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(openxlsx)
library(tidyverse)
library(RColorBrewer)

opt <- NULL
opt$Matrix <- "116EVssample.counts.txt"
opt$Anno <- "gene_annotation.xlsx"
opt$Group <- "116samples分组情况.xls"
opt$Prefix <- "Samples_RNA_distribution"
opt$outdir <- "."

command <- matrix(
    c("Matrix","D",1,"character",
      "Anno","A",1,"character",
      "Group","G",1,"character",
      "Prefix","P",1,"character",
      "outdir","O",2,"character")
      ,byrow=T,ncol=4
)
opt <- getopt(command)

if(!dir.exists(opt$outdir)){
    dir.create(opt$outdir)
}

da <- read.table(opt$Matrix,sep="\t",header=T,stringsAsFactors=F,check.names=F)
anno <- read.xlsx(opt$Anno)
groups <- read.table(opt$Group,header=T,sep="\t",stringsAsFactors=F,check.names=F)

info <- merge(da,anno,by.x="Geneid",by.y="ensembl_gene_id",all.x=T)

IGtype <- c("IG_C_gene","IG_D_gene","IG_J_gene","IG_V_gene")
lncRNAtype <- c("3prime_overlapping_ncRNA","antisense_RNA","macro_lncRNA","non_coding","bidirectional_promoter_lncRNA","sense_intronic","sense_overlapping","lincRNA")
ncRNAtype <- c("miRNA","misc_RNA","piRNA","rRNA","siRNA","snRNA","snoRNA","tRNA","vaultRNA","Mt_rRNA","Mt_tRNA")
proteinCodingtype <- c("protein_coding")
Pseudogenetype <- c("IG_pseudogene","IG_C_pseudogene","IG_J_pseudogene","IG_V_pseudogene","polymorphic_pseudogene","processed_pseudogene","pseudogene","transcribed_processed_pseudogene","transcribed_unitary_pseudogene","transcribed_unprocessed_pseudogene","translated_processed_pseudogene","TR_J_pseudogene","TR_V_pseudogene","unitary_pseudogene","unprocessed_pseudogene")
TRtype <- c("TR_C_gene","TR_D_gene","TR_J_gene","TR_V_gene")

typeList <- list(IG_gene=IGtype,lncRNA=lncRNAtype,ncRNA=ncRNAtype,proteinCoding=proteinCodingtype,pseudogene=Pseudogenetype,TR_gene=TRtype)

get_type_count <- function(x,da){
	types <- da$gene_biotype[which(x!=0)]
	typeCount <- sapply(typeList,function(x){return(length(which(types %in% x)))})
	others <- length(types)-sum(typeCount)
	names(others) <- "others"
	typeCount <- c(typeCount,others)
	return(typeCount)
}

type_count <- do.call(rbind,lapply(info[7:(ncol(info)-3)],function(x) get_type_count(x,info)))

type_count_anno <- type_count %>% as.data.frame() %>% rownames_to_column() %>% rename(SampleID=rowname) %>% merge(groups,by.x="SampleID",by.y="SampleID",all.y=T)

type_count_anno_long <- melt(type_count_anno,id.vars=colnames(type_count_anno)[-(1:(ncol(type_count_anno)-4))],measure.vars = colnames(type_count_anno)[2:(ncol(type_count_anno)-4)],variable.name="RNAtype",value.name="Counts")

write.table(type_count_anno,file=paste0(opt$outdir,"/",opt$Prefix,'.txt'),sep="\t",quote=F,row.names=F,col.names=T)

Boxplot <- function(type_count_anno_long,cancer,mode="All",outdir){
	case <- filter(type_count_anno_long,grepl(paste0('^',cancer,'.*','Cancer$'),Group2)) %>% mutate(condition=paste0(cancer,'_Case'))
	if(mode=='All'){
                ctrl <- filter(type_count_anno_long,grepl(paste0('^',cancer,'.*','control$'),Group2) | grepl('^Shared.*control$',Group2)) %>% mutate(condition=paste0(cancer,"_All_Ctrl"))
        }else if(mode=='CtrlOnly'){
                ctrl <- filter(type_count_anno_long,grepl(paste0('^',cancer,'.*','control$'),Group2)) %>% mutate(condition=paste0(cancer,'_Ctrl'))
        }else if(mode=='ShareOnly'){
                ctrl <- filter(type_count_anno_long,grepl('^Shared.*control$',Group2)) %>% mutate(condition='SharedCtrl')
        }else{
                stop("Wrong mode")
        }
	type_count <- as.data.frame(rbind(case,ctrl))
	p <- ggboxplot(type_count,x="RNAtype",y="Counts",fill="condition",border=NULL,font.label=list(size=5)) + stat_compare_means(aes(group = condition), label = "p.format") + theme(text = element_text(size=10),axis.text.x = element_text(angle=45, hjust=1))
	ggsave(p,file=paste0(outdir,'/',cancer,'_sample_RNA_types_boxplot.pdf'))
}

cancerTypes <- c("Bladder","Kidney")
sapply(cancerTypes,function(x) Boxplot(type_count_anno_long,x,mode="All",opt$outdir))
