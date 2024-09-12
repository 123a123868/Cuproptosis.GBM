# 补充图片：
# 免疫细胞直接的相互作用网络图

# 免疫细胞的相关性网络图

# 设置
#设置
rm(list=ls())
options(stringsAsFactors = F)
library(reshape2)
library(corrplot)
library(plyr)
library(igraph) #用于绘制网络图
# 颜色
poscol <- "#FB9A99" #正相关用红色连线
negcol <- "#C6DBEF" #负相关用蓝色连线
mycol <- c("#FDBF6F", "#1F78B4", "#E31A1C", "#8C510A") #cluster的颜色，如果有更多类，就给更多的颜色

# 免疫细胞矩阵
# 免疫细胞矩阵
input <- "Figure4/out_data/ssgseaOut.txt"  
immune <- read.table( input ,sep="\t",header=T,row.names=1,check.names=F)
immune[1:4, 1:4]
input_data = as.data.frame(t(immune))
input_data[1:3,1:3]


# 计算相关性，然后分了4群
# 这个4群可以自己设置
# options(scipen = 100,digits = 4)
corr <- cor(input_data, method = "spearman")
corrplot(corr,title = "", 
         method = "pie", #或"circle" (default), "square", "ellipse", "number", "pie", "shade" and "color"
         outline = T, addgrid.col = "darkgray", 
         order="hclust", addrect = 4, #hclust聚为4类，根据数据的具体情况调整
         mar = c(4,0,4,0), #撑大画布，让细胞名显示完全
         rect.col = "black", rect.lwd = 5, cl.pos = "b", 
         tl.col = "black", tl.cex = 1.08, cl.cex = 1.5, tl.srt=60)

cor.mtest <- function(corr, ...) {
  corr <- as.matrix(corr)
  n <- ncol(corr)
  p.corr <- matrix(NA, n, n)
  diag(p.corr) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(corr[, i],method = "spearman", corr[, j], ...)
      p.corr[i, j] <- p.corr[j, i] <- tmp$p.value
    }
  }
  colnames(p.corr) <- rownames(p.corr) <- colnames(corr)
  p.corr
}
p.corr <- cor.mtest(input_data)  # 计算p值
head(p.corr[, 1:5])


#合并相关系数和P值
rr <- as.data.frame(corr);
rr$ID <- rownames(rr)
cor <- melt(rr,"ID",value.name = "cor"); #head(cor)
pp <- as.data.frame(p.corr);
pp$ID <- rownames(pp)
pvalue <- melt(pp,"ID",value.name = "pvalue"); #head(pvalue)
colnames(pvalue) <- c("from","to","pvalue")
corpvlue <- cbind(pvalue, cor)
head(corpvlue)
corpvlue <- corpvlue[, -c(4:5)]
head(corpvlue)
dim(corpvlue)


#去掉相关性较弱的连接
corpvlue <- corpvlue[corpvlue$pvalue < 0.01,] #只保留pvalue < 0.01的
dim(corpvlue)
# 计算权重，线条的粗细
corpvlue$weight <- corpvlue$pvalue
corpvlue$weight <- -log10(corpvlue$weight)
head(corpvlue)
#去掉相关系数为1，也就是两个相同变量之间的连接
corpvlue <- corpvlue[!corpvlue$cor==1,]
dim(corpvlue)
#去掉相关系数一样的连接--也就是重复计算的连接
summary(duplicated(corpvlue$weight))
corpvlue <- corpvlue[!duplicated(corpvlue$weight),]
dim(corpvlue)
#相关系数的正负用不同颜色表示
corpvlue$color <- ifelse(corpvlue$cor<0, negcol, poscol)
#保存到文件，便于查看
write.csv(corpvlue, "10-F4/output_links.csv")

# 分成4类
cellcluster <- as.data.frame(t(input_data))
#cellcluster[1:5,1:5]
hc <- hclust(dist((cellcluster)))
hcd <- as.dendrogram(hc)
(clus4 <- cutree(hc, 4)) #分4类
A <- as.character(rownames(as.data.frame(subset(clus4,clus4==1))))
B <- as.character(rownames(as.data.frame(subset(clus4,clus4==2))))
C <- as.character(rownames(as.data.frame(subset(clus4,clus4==3))))
D <- as.character(rownames(as.data.frame(subset(clus4,clus4==4))))
cls <- list(A,B,C,D)

nodes <- as.data.frame(unlist(cls))
# 这里得自己修改，看看每群有几个--cls
nodes$type <- c(rep("A",5),rep("B",5),rep("C",5),rep("D",9))
names(nodes) <- c("media","type.label")

#以hclust的结果为基础，调整部分细胞所属的cluster
# nodes$type.label[nodes$media=="T cells follicular helper"] <- "B"
# nodes$type.label[nodes$media=="B cells naive"] <- "A"
# nodes$type.label[nodes$media=="T cells CD4 naive"] <- "A"
# nodes$type.label[nodes$media=="Plasma cells"] <- "A"
# nodes$type.label[nodes$media=="Dendritic cells resting"] <- "C"
# nodes$type.label[nodes$media=="Eosinophils"] <- "C"
# nodes$type.label[nodes$media=="Mast cells resting"] <- "A"
nodes <- as.data.frame(nodes)
nodes$media <- as.character(nodes$media)
nodes

# 合并生存分析的数据和细胞分类的数据
# summary(nodes$media %in% bb$ID) #检查细胞名是否一致
# nodes <- merge(nodes, bb, by.x = "media", "ID", all.x = T, all.y = T) #按细胞名merge
# nodes$Fraction <- abs(nodes$weight_HR)
nodes$id <- paste("S", 01:24, sep = "")
nodes <- nodes[order(nodes$type.label),]
nodes <- nodes[,c(ncol(nodes),1:ncol(nodes)-1)]
nodes <- nodes[order(nodes$type.label),]
nodes


#建立nodes和links的连接id，把细胞名换成ID
paste0("'",nodes$media,"'","=","'",nodes$id,"'",collapse = ",")
corpvlue$from <- revalue(corpvlue$from,
                         c('B cells'='S1','Cytotoxic cells'='S2','aDC'='S3','Mast cells'='S4','Neutrophils'='S5',
                           'T cells'='S6','Th17 cells'='S7','TReg'='S8','DC'='S9','pDC'='S10','T helper cells'='S11',
                           'CD8 T cells'='S12','Tgd'='S13','NK CD56dim cells'='S14','Macrophages'='S15','Tcm'='S16',
                           'Tem'='S17','Th1 cells'='S18','Th2 cells'='S19','TFH'='S20','NK cells'='S21',
                           'NK CD56bright cells'='S22','iDC'='S23','Eosinophils'='S24'))


corpvlue$to <- revalue(corpvlue$to,c('B cells'='S1','Cytotoxic cells'='S2','aDC'='S3','Mast cells'='S4','Neutrophils'='S5',
                                     'T cells'='S6','Th17 cells'='S7','TReg'='S8','DC'='S9','pDC'='S10','T helper cells'='S11',
                                     'CD8 T cells'='S12','Tgd'='S13','NK CD56dim cells'='S14','Macrophages'='S15','Tcm'='S16',
                                     'Tem'='S17','Th1 cells'='S18','Th2 cells'='S19','TFH'='S20','NK cells'='S21',
                                     'NK CD56bright cells'='S22','iDC'='S23','Eosinophils'='S24'))
(links <- corpvlue)


#利用nodes和links构建网络的input文件
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 



#####################################################################################
# Generate colors based on cell clusters:
V(net)$color <- revalue(nodes$type.label,c("A"=mycol[1],"B"=mycol[2],"C"=mycol[3],"D"=mycol[4]))
# Compute node degrees (#links) and use that to set node size:
# Set edge width based on weight-log10(p_value):
V(net)$size <- (1 + V(net)$weight)*3 #节点圆的大小，可根据自己的数据再做调整
V(net)$label <- V(net)$media #设置标签
E(net)$arrow.mode <- 0 #不需要箭头
E(net)$edge.color <- "tomato" # tomato gray80
E(net)$width <- 1+E(net)$weight/6  #连接之间权重

pdf("Figure4/out_polt//Immune_network.pdf", width = 9.75, height = 8.78 )
plot(net,
     layout=layout_in_circle, #按圆圈布局
     edge.curved=.2, #画弯曲的连线
     vertex.label.color=V(net)$color, #细胞名的颜色
     vertex.label.dist= -2, #标签和节点的位置错开，后期还是要用AI调整
     edge.color=links$color)

#cluster的图例
legend("topright", #图例的位置
       c("Cell cluster-A", "Cell cluster-B", "Cell cluster-C", "Cell cluster-D"),
       pch=21, col="black", pt.bg=mycol, pt.cex=3,
       cex=1.3, bty="n", ncol=1)

#节点圆大小的图例，参考了FigureYa75base_volcano
f <- c(0.05, 0.001, 0.00001, 0.00000001)
s <- sqrt(abs(log10(f)))*3
legend("bottomright", 
       inset=c(0,-.1), #向下移
       legend=f, text.width = .2, 
       title = "logrank test, P value", title.adj = -.5,
       pch=21, pt.cex=s, bty='n',
       horiz = TRUE, #横向排列
       col = "black")

#连线的图例
legend("bottomright",
       c("Positive correlation with P < 0.0001", 
         "Negative correlation with P < 0.0001"),
       col = c(poscol, negcol), bty="n", 
       cex = 1, lty = 1, lwd = 5)

dev.off()


#######################################################################
# 计算8个基因与图A中22群免疫细胞之间的相关性

# 设置
#设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)
library(tinyarray)
library(tidyverse)

# 的表达矩阵
# 读入数据
# 读入数据
rm(list = ls()) 
# 读入10铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene

# 导入表达矩阵数据
# 读入数据
load("Figure3/out_data/expr.gbm.Rdata")
data_gene = expr.gbm[cu_gene,]
data_gene = as.data.frame(t(data_gene))

# 免疫细胞矩阵
input <- "Figure4/out_data/ssgseaOut.txt"  
immune <- read.table( input ,sep="\t",header=T,row.names=1,check.names=F)
immune[1:4, 1:4]
input_data = as.data.frame(t(immune))
input_data[1:3,1:3]
identical(rownames(data_gene),rownames(input_data))

# 计算相关性，并且画热图
sig_gene <- cu_gene
library(psych)
x <- data_gene
y <- input_data
library(psych)
d <- corr.test(x,y,use="complete",method = 'spearman')
r <- d$r
p <- d$p
library(ggcorrplot)
ggcorrplot(t(d$r), show.legend = T, colors = c("#6D9EC1", "white", "#E46726"),
           p.mat = t(d$p.adj), digits = 2,  sig.level = 0.05,insig = 'blank',lab = F
)
ggsave("Figure4/out_polt//cor-Cibersort-immune.pdf")
###############










