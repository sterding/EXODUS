library(tidyverse)
library(DESeq2)
library(edgeR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(BiocParallel)
register(MulticoreParam(8))

Info <- read.table("116EVssample.counts.txt",header=T,row.names=1,sep="\t",stringsAsFactors=F,check.names=F)
countMatrix <- Info[,6:ncol(Info)]
groups <- read.table("116samples分组情况.xls",header=T,sep="\t",stringsAsFactors=F)


DE_analysis <- function(mat,groups,cancer,mode='All',frac=0.8,seed,foldChangeThreshold=0.585,method="DEseq2",pathFun=pathwary_analysis){
	## group samples by cancer type and mode
	case <- filter(groups,grepl(paste0('^',cancer,'.*','Cancer$'),Group2)) %>% mutate(condition=paste0(cancer,'_Case'))
	if(mode=='All'){
		ctrl <- filter(groups,grepl(paste0('^',cancer,'.*','control$'),Group2) | grepl('^Shared.*control$',Group2)) %>% mutate(condition=paste0(cancer,"_AllCtrl"))
	}else if(mode=='CtrlOnly'){
		ctrl <- filter(groups,grepl(paste0('^',cancer,'.*','control$'),Group2)) %>% mutate(condition=paste0(cancer,'_Ctrl'))
	}else if(mode=='ShareOnly'){
		ctrl <- filter(groups,grepl('^Shared.*control$',Group2)) %>% mutate(condition='SharedCtrl')
	}else{
		stop("Wrong mode")
	}

	## stratified sampling for train and test data set
	set.seed(seed)
	trainCase <- sample_frac(case,frac)
	trainCtrl <- sample_frac(ctrl,frac)
	testCase <- filter(case,!case$SampleID %in% trainCase$SampleID)
	testCtrl <- filter(ctrl,!ctrl$SampleID %in% trainCtrl$SampleID)
	trainColdata <- rbind(trainCase,trainCtrl) %>% dplyr::select(SampleID,condition) %>% mutate_at("condition",factor) %>% column_to_rownames("SampleID")
	testColdata <- rbind(testCase,testCtrl) %>% dplyr::select(SampleID,condition) %>% mutate_at("condition",factor) %>% column_to_rownames("SampleID")
	countdata <- mat %>% dplyr::select(rownames(trainColdata))
	group <- levels(trainColdata$condition)

	## DE analysis
	if(method=='DEseq2'){
		dds <- DESeqDataSetFromMatrix(countData=countdata,colData=trainColdata,design=~condition)
		dds <- DESeq(dds)
		res <- results(dds,contrast=c("condition",group[grep('Case$',group,perl=T)],group[grep('Ctrl$',group,perl=T)])) %>% as.data.frame %>% rownames_to_column() %>% rename(c('rowname'='EnsemblID')) %>% arrange(pvalue)
		resFilterPadj <- filter(res,padj < 0.05 & abs(log2FoldChange)>foldChangeThreshold) %>% arrange(padj)
		resFilterPvalue <- filter(res,pvalue < 0.05 & abs(log2FoldChange)>foldChangeThreshold) %>% arrange(pvalue)
	}else if(method=='edgeR_quasi-likelihood'){
		Treat <- factor(trainColdata$condition,levels=(c(group[grep('Ctrl$',group,perl=T)],group[grep('Case$',group,perl=T)])))
		DGE <-  DGEList(counts=countdata, group=Treat)
		#KEEP <-  filterByExpr(DGE)
		KEEP <-  rowSums(cpm(DGE)>1)>=3
		DGE <- DGE[KEEP, , keep.lib.sizes=FALSE]
		DGEnorm <- calcNormFactors(DGE)
		design <- model.matrix(~Treat)
		DGEdisp <- estimateDisp(DGEnorm, design, robust=TRUE)
		fit <- glmQLFit(DGEdisp, design, robust=TRUE)
		qlf <- glmQLFTest(fit)
		res <- qlf$table
		res$FDR <- p.adjust(res$PValue, method="BH")
		resFilterPadj <- res %>% rownames_to_column() %>% filter(FDR < 0.05 & abs(logFC)>foldChangeThreshold) %>% rename(c('logFC'='log2FoldChange','rowname'='EnsemblID')) %>% arrange(FDR)
		resFilterPvalue <- res %>% rownames_to_column() %>% filter(PValue < 0.05 & abs(logFC)>foldChangeThreshold) %>% rename(c('logFC'='log2FoldChange','rowname'='EnsemblID')) %>% arrange(PValue)
	}else if(method=='edgeR_likelihood'){
		Treat <- factor(trainColdata$condition,levels=(c(group[grep('Ctrl$',group,perl=T)],group[grep('Case$',group,perl=T)])))
		DGE <-  DGEList(counts=countdata, group=Treat)
		KEEP <-  rowSums(cpm(DGE)>1)>=3
		DGE <- DGE[KEEP, , keep.lib.sizes=FALSE]
		DGEnorm <- calcNormFactors(DGE)
		design <- model.matrix(~Treat)
		DGEdisp <- estimateDisp(DGEnorm, design, robust=TRUE)
		fit <- glmFit(DGEdisp,design)
		lrt <- glmLRT(fit)
		res <- lrt$table
		res$FDR <- p.adjust(res$PValue, method="BH")
		resFilterPadj <- res %>% rownames_to_column() %>% filter(FDR < 0.05 & abs(logFC)>foldChangeThreshold) %>% rename(c('logFC'='log2FoldChange','rowname'='EnsemblID')) %>% arrange(FDR)
		resFilterPvalue <- res %>% rownames_to_column() %>% filter(PValue < 0.05 & abs(logFC)>foldChangeThreshold) %>% rename(c('logFC'='log2FoldChange','rowname'='EnsemblID')) %>% arrange(PValue)
	}else if(method=='edgeR_classic'){
		Treat <- factor(trainColdata$condition,levels=(c(group[grep('Ctrl$',group,perl=T)],group[grep('Case$',group,perl=T)])))
		DGE <-  DGEList(counts=countdata, group=Treat)
		KEEP <-  rowSums(cpm(DGE)>1)>=3
		DGE <- DGE[KEEP, , keep.lib.sizes=FALSE]
		DGEnorm <- calcNormFactors(DGE)
		DGEdisp <- estimateDisp(DGEnorm)
		et <- exactTest(DGEdisp)
		res <- et$table
		res$FDR <- p.adjust(res$PValue, method="BH")
		resFilterPadj <- res %>% rownames_to_column() %>% filter(FDR < 0.05 & abs(logFC)>foldChangeThreshold) %>% rename(c('logFC'='log2FoldChange','rowname'='EnsemblID')) %>% arrange(FDR)
		resFilterPvalue <- res %>% rownames_to_column() %>% filter(PValue < 0.05 & abs(logFC)>foldChangeThreshold) %>% rename(c('logFC'='log2FoldChange','rowname'='EnsemblID')) %>% arrange(PValue)
	}else{
		stop("Wrong DE method")
	}
	trainDEgeneMatrixPadj <- mat %>% rownames_to_column() %>% rename(c('rowname'='EnsemblID')) %>% filter(EnsemblID %in% resFilterPadj$EnsemblID) %>% dplyr::select(EnsemblID,rownames(trainColdata))
	testDEgeneMatrixPadj <- mat %>% rownames_to_column() %>% rename(c('rowname'='EnsemblID')) %>% filter(EnsemblID %in% resFilterPadj$EnsemblID) %>% dplyr::select(EnsemblID,rownames(testColdata))
	trainDEgeneMatrixPval <- mat %>% rownames_to_column() %>% rename(c('rowname'='EnsemblID')) %>% filter(EnsemblID %in% resFilterPvalue$EnsemblID) %>% dplyr::select(EnsemblID,rownames(trainColdata))
	testDEgeneMatrixPval <- mat %>% rownames_to_column() %>% rename(c('rowname'='EnsemblID')) %>% filter(EnsemblID %in% resFilterPvalue$EnsemblID) %>% dplyr::select(EnsemblID,rownames(testColdata))
	prefix <- paste0(group[grep('Case$',group,perl=T)],'_VS_',group[grep('Ctrl$',group,perl=T)])
	dir <- paste0(prefix,'_',method)
	dir.create(dir)
	write.table(rbind(trainCase,trainCtrl),file=paste0(dir,'/',prefix,'_train_dataSet_group.tsv'),sep="\t",quote=F,col.names=T,row.names=F)
	write.table(rbind(testCase,testCtrl),file=paste0(dir,'/',prefix,'_test_dataSet_group.tsv'),sep="\t",quote=F,col.names=T,row.names=F)
	write.table(res,file=paste0(dir,'/',prefix,'_',method,'_DE_all.tsv'),sep="\t",quote=F,col.names=T,row.names=F)
	write.table(resFilterPadj,file=paste0(dir,'/',prefix,'_',method,'_DE_log2foldChange',foldChangeThreshold,'Padj0.05.tsv'),sep="\t",quote=F,col.names=T,row.names=F)
	write.table(resFilterPvalue,file=paste0(dir,'/',prefix,'_',method,'_DE_log2foldChange',foldChangeThreshold,'Pvalue0.05.tsv'),sep="\t",quote=F,col.names=T,row.names=F)
	write.table(trainDEgeneMatrixPadj,file=paste0(dir,'/',prefix,'_',method,'_trainDEgeneMatrix_padj.tsv'),sep="\t",quote=F,col.names=T,row.names=F)
	write.table(testDEgeneMatrixPadj,file=paste0(dir,'/',prefix,'_',method,'_testDEgeneMatrix_padj.tsv'),sep="\t",quote=F,col.names=T,row.names=F)
	write.table(trainDEgeneMatrixPval,file=paste0(dir,'/',prefix,'_',method,'_trainDEgeneMatrix_pval.tsv'),sep="\t",quote=F,col.names=T,row.names=F)
	write.table(testDEgeneMatrixPval,file=paste0(dir,'/',prefix,'_',method,'_testDEgeneMatrix_pval.tsv'),sep="\t",quote=F,col.names=T,row.names=F)

	## pathway analysis
	regulateStyles <- c('All','Up','Down')
	sapply(regulateStyles,function(x){pathFun(resFilterPvalue,x,prefix,dir)})

	## save Rdata
	save(list=ls(),file=paste0(dir,'/',prefix,'_',method,'.Rdata'))
}

pathwary_analysis <- function(DEres,regulateStyle='All',prefix,dir){
	if(regulateStyle=='All'){
		DEres <- DEres %>% arrange(desc(log2FoldChange))
	}else if(regulateStyle=='Up'){
		DEres <- DEres %>% filter(log2FoldChange>0) %>% arrange(desc(log2FoldChange))
	}else if(regulateStyle=='Down'){
		DEres <- DEres %>% filter(log2FoldChange<0) %>% arrange(desc(log2FoldChange))
	}else{
		stop("Wrong style")
	}
	pathDir <- paste0(dir,'/pathway_result_with_regulate_style_',regulateStyle)
	dir.create(pathDir)
	map <- bitr(DEres$EnsemblID,fromType="ENSEMBL",toType=c("ENTREZID","SYMBOL"),OrgDb="org.Hs.eg.db")
	mapTable <- data.frame(cbind(map,DEres[match(map$ENSEMBL,DEres$EnsemblID),]$log2FoldChange))
	colnames(mapTable)[4] <- 'log2FoldChange'
	mapTable <- mapTable[!duplicated(mapTable$ENTREZID),]
	geneList <- mapTable$log2FoldChange %>% set_names(mapTable$ENTREZID)
	#test.gene <- data.frame(ensembl = DE$EnsemblID, entrez = NA, score = DE$log2FoldChange, stringsAsFactors = F)      			
	#test.gene$ensembl <- substr(test.gene$ensembl,1,15)
	#test.gene$entrez <- mapIds(org.Hs.eg.db, keys=test.gene$ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
	#test.gene <- test.gene[!is.na(test.gene$entrez),]
	ego <- enrichGO(gene=mapTable$ENTREZID,OrgDb=org.Hs.eg.db,ont="ALL",pAdjustMethod = "BH",pvalueCutoff=0.05)
	ekk <- enrichKEGG(gene=mapTable$ENTREZID,organism="hsa",pvalueCutoff=0.05)
	gsego <- gseGO(geneList = geneList,
          OrgDb        = org.Hs.eg.db,
          ont          = "All",
          nPerm        = 1000,
          minGSSize    = 100,
          maxGSSize    = 500,
          pvalueCutoff = 0.05,
          verbose      = FALSE)
	tryCatch({gsekk <- gseKEGG(geneList = geneList,
           organism     = 'hsa',
           nPerm        = 1000,
           minGSSize    = 120,
           pvalueCutoff = 0.05,
           verbose      = FALSE)},
           error = function(e){cat("ERROR :",conditionMessage(e),"\n")}
    )
    if(dim(ego)[1] != 0){
	write.table(ego,paste0(pathDir,'/',prefix,'_go_ora.tsv'),sep="\t",quote=F,col.names=T,row.names=F)
	library(enrichplot)
	#pdf(paste0(pathDir,'/',prefix,'_go_ora.pdf'),width=8,height=6)
	p1 <- barplot(ego,font.size=10,showCategory=10)
	#dev.off()
	ggsave(p1,file=paste0(pathDir,'/',prefix,'_go_ora.pdf'),width=8,height=6)
    }
    if(dim(ekk)[1] != 0){
	write.table(ekk,paste0(pathDir,'/',prefix,'_kegg_ora.tsv'),sep="\t",quote=F,col.names=T,row.names=F)
	library(enrichplot)
	#pdf(paste0(pathDir,'/',prefix,'_kegg_ora.pdf'),width=8,height=6)
	p2 <- barplot(ekk,font.size=10,showCategory=10)
	#dev.off()
	ggsave(p2,file=paste0(pathDir,'/',prefix,'_kegg_ora.pdf'),width=8,height=6)
    }
    if(dim(gsego)[1] != 0){
	write.table(gsego,paste0(pathDir,'/',prefix,'_go_gse.tsv'),sep="\t",quote=F,col.names=T,row.names=F)
	library(enrichplot)
	#pdf(paste0(pathDir,'/',prefix,'_go_gse.pdf'),width=8,height=6)
	p3 <- dotplot(gsego,font.size=10,showCategory=10)
	#dev.off()
	ggsave(p3,file=paste0(pathDir,'/',prefix,'_go_gse.pdf'),width=8,height=6)
    }
    if(exists('gsekk') & dim(gsekk)[1] !=0){
	write.table(gsekk,paste0(pathDir,'/',prefix,'_kegg_gse.tsv'),sep="\t",quote=F,col.names=T,row.names=F)
	library(enrichplot)
	#pdf(paste0(pathDir,'/',prefix,'_kegg_gse.pdf'),width=8,height=6)
	p4 <- dotplot(gsekk,font.size=10,showCategory=10)
	#dev.off()
	ggsave(p4,file=paste0(pathDir,'/',prefix,'_kegg_gse.pdf'),width=8,height=6)
    }
}

Cancers <- c("Bladder","Kidney")
Methods <- c("DEseq2","edgeR_quasi-likelihood","edgeR_likelihood","edgeR_classic")

for(i in 1:length(Cancers)){
	sapply(Methods,function(x){DE_analysis(mat=countMatrix,groups=groups,cancer=Cancers[i],mode='All',frac=0.8,seed=i,foldChangeThreshold=0.585,method=x,pathFun=pathwary_analysis)})
}
