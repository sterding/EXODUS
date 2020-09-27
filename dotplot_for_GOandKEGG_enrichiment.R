#######################
library(ggplot2)
#######################

#########KEGG#########
#0 usage
data <- "KEGG.txt"

#1 data process
data <- read.table(data, head=T, sep="\t", check.names=F)
dat <- data[data$rawP < 0.1, ]
dat <- dat[order(-dat$FoldChange, dat$rawP), ] # ½µÐò + ÉýÐò ÅÅÁÐ
dat$name <- as.character(dat$Term_Name)
for (i in 1:nrow(dat)){
  if(nchar(dat$name[i]) >= 60){
    dat$name[i] <- paste(substring(dat$name[i], 1, 60), "..")
  }
}
dat$name <- factor(dat$name, levels=dat$name)

#2 theme
mytheme <- theme_bw() + 
  theme(text=element_text(colour="black"),
        plot.title=element_text(size=rel(1.5),lineheight=0.7,hjust=0.5),
        axis.text=element_text(size=rel(1)),
        axis.title=element_text(size=rel(1.4),lineheight=0.4,hjust=0.5),
        panel.border=element_rect(size=0.3,colour="black"),
        axis.line=element_blank(),
        axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),
)

#3 plot
p<-ggplot(dat, aes(x=FoldChange, y=name, size=Observed, colour=rawP)) + 
   geom_point() + 
   mytheme + ylab(NULL) +
   scale_color_continuous(high="#000099", low="#FF0000", space='Lab', name="significance level") +
   guides(size=guide_legend(title="observed genes"))
p + ggtitle(paste("\nthe most enriched pathways\n"))

##########GO#########
input <- "GO.txt"
data <- read.table(input, head=T, sep="\t", check.names=F)
dat <- data[data$rawP < 0.05 , ]
dat <- dat[order(dat$class, -dat$FoldChange, dat$rawP), ] # ½µÐò + ÉýÐò ÅÅÁÐ
dat$name <- as.character(dat$Term_Name)
for (i in 1:nrow(dat)){
  if(nchar(dat$name[i]) >= 60){
    dat$name[i] <- paste(substring(dat$name[i], 1, 60), "..", sep="")
  }
}
dat$name <- factor(dat$name, levels=dat$name)
out <- dat
######

#2 theme
mytheme<-theme_bw() + 
  theme(text=element_text(colour="black"),
        plot.title=element_text(size=rel(1.5),lineheight=0.7,hjust=0.5),
        axis.text=element_text(size=rel(1)),
        axis.title=element_text(size=rel(1.4),lineheight=0.4,hjust=0.5),
        panel.border=element_rect(size=0.3,colour="black"),
        axis.line=element_blank(),
        axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),
)

#3 plot
dat <- out
p <- ggplot(dat, aes(x=FoldChange, y=name, size=Observed, colour=rawP)) + 
     geom_point() + 
     mytheme + ylab(NULL) +
     scale_color_continuous(high="#000099", low="#FF0000", space='Lab', name="significance level") +
     guides(size=guide_legend(title="Observed genes")) +
     ggtitle(paste("\nthe most enriched GO terms\n")) + 
     facet_grid(class~., scale="free_y", shrink = F, space="free_y" )  # ·ÖÀ¸
p
##############END