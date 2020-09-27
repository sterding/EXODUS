# Rscript in.R input background_total gene_total Adjustment out
args    <- commandArgs(TRUE)
infile  <- args[1]
N 	<- as.numeric(args[2])
T	<- as.numeric(args[3])
adj     <- args[4]
outfile <- args[5]

rt  <- read.table(infile, sep="\t", head=T)

dat <- rt[ rt[,3]>=1, ]

mat <- dat[,2:3]
mat <- cbind(mat,0)
mat <- cbind(mat,0)

for(i in 1:nrow(mat)){
	m <- matrix( c(N,T,mat[i,1],mat[i,2]), nrow=2, dimnames=list(c("Background","Input"), c("Total","Term")) )
	mat[i,3] <- fisher.test(m, alternative="two.sided")$p.value
}

mat[,4] <- p.adjust( mat[,3], method=adj, n=length(mat[,3]) )

out <- cbind(dat,mat)[, c(1:5,9,10,6)]
colnames(out) <- c("Term","Category","Observed","Expected","FoldChange","rawP","adjP","Genelist")

out <- out[order(out[,7]),]

write.table(out, file=outfile, quote=F, sep = "\t", row.names = F)
