#先获取矩阵数据

#目的：IMvigor210CoreBiologies  免疫治疗的分析
#用的包：IMvigor210CoreBiologies 

#设置
rm(list = ls()) 
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
library(IOBR) 
library(IMvigor210CoreBiologies)

# 读入数据
load("Figure7/out_data//IMvigor_count.Rdata")
dim(expr)
expr[1:4,1:4]
head(pdata)


# cc分型
# 读入10铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
expr.cu = expr[cu_gene,]

#sweep函数减去中位数进行标准化
d = sweep(expr.cu,1, apply(expr.cu,1,median,na.rm=T))
d = as.matrix(expr.cu)


#2、CC分型
setwd("Figure7/out_data/")
dir.create('IM_CC-result')
title="./IM_CC-result"
library(ConsensusClusterPlus)
res <- ConsensusClusterPlus(d, maxK = 7, reps = 1000, pItem = 0.8, 
                            pFeature = 1,corUse = "complete.obs", seed=123456, plot="pdf", writeTable=T,
                            title=title)
setwd("../../")
getwd()
######################################################
# 可视化

# 设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)
library(data.table)
library(dplyr)
library(futile.logger) 
library(glmSparseNet)
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)


# 读入数据
load("Figure7/out_data//IMvigor_count.Rdata")
dim(expr)
expr[1:4,1:4]
head(pdata)
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
expr.cu = expr[cu_gene,]

# 读入CC分型信息
pd.CC = fread( 'Figure7/out_data/IM_CC-result/IM_CC-result.k=2.consensusClass.csv',
               header = F, data.table = F)
colnames(pd.CC) = c('sample','cluster')
pd.CC$cluster = paste0('C',pd.CC$cluster) 

#匹配合并一下：
rownames(pdata) = pdata$ID
pdata = pdata[pd.CC$sample,]
identical(rownames(pdata),pd.CC$sample)
pdata$CC_cluter = pd.CC$cluster

# 绘图数据
library(dplyr)
library(plyr)
a <- data.frame(table(pdata$CC_cluter,pdata$BOR_binary))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
a$Freq

pvalue <- chisq.test(c(49,185,19,45,ncol=2))$p.value #卡方检验
library(plyr)
ggplot(a,aes(Var1,percent,fill=Var2))+
  geom_bar(stat="identity",position = position_stack())+
  scale_fill_manual(values = c("#008ECA","#DB423E"),label=c("CR/PR","SD/PD"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="type",y="Percent Weidght",
       fill="")+
  geom_text(aes(label=label),vjust=3,size=6,color="black")+
  annotate(geom = "text",
           cex=6,
           x=1.5, y=105, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
ggsave(filename = "Figure7/out_polt/F7_C.pdf")

