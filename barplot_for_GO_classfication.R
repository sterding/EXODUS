### read data
data <- read.table("bladder.go_class.txt", sep="\t", header=T)
data <- read.table("kidney.go_class.txt", sep="\t", header=T)
head(data)

### handle it, 将term name设定为factor，即可按照顺序输出
term_order = factor(as.integer(rownames(data)), labels=data$name)
term_colos = c("#66C3A5", "#8DA1CB", "#FD8D62")

### plot 
library(ggplot2)
ggplot(data=data, 
  aes(x=term_order, y=gene, fill=class)) +
  geom_bar(stat="identity", width=0.8)  + 
  scale_fill_manual(values= term_colos) + 
  theme_bw() +
  xlab("GO term") + 
  ylab("Num of Genes") + 
  labs(title = "goslim summary") + 
  theme(axis.text.x=element_text(face = "bold", color="gray50", angle=70, vjust=1, hjust=1)
) 
###
