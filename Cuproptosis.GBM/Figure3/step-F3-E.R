# C1和C2的 10CU基因的热图表达

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)

#####################
#1、处理数据
#    CGGA数据库
#####################
# 表达矩阵
library(data.table)
a = fread( 'Input_data/CGGA_GBM/CGGA.mRNAseq_693.RSEM-genes.20200506.txt/CGGA.mRNAseq_693.RSEM-genes.20200506.txt',
           header = T, data.table = F)
dim(a) 
a[1:4,1:4]
rownames(a) = a$Gene_Name
expr.a = a[,-1]
expr.a[1:4,1:4]
boxplot(expr.a[,1:8],las=2)
expr.a = log2(expr.a+1)


#读入CC分型信息
pd.CC = fread( 'Figure3/CC-result/CC-result.k=2.consensusClass.csv',
               header = F, data.table = F)
colnames(pd.CC) = c('sample','cluster')
pd.CC$cluster = paste0('C',pd.CC$cluster) 
# 匹配取交集
expr.gbm = expr.a[,pd.CC$sample]
dim(expr.gbm)
identical(colnames(expr.gbm),pd.CC$sample)
# 保存数据
save(expr.gbm,file = "Figure3/out_data/expr.gbm.Rdata")
save(pd.CC,file = "Figure3/out_data/pd.CC.Rdata")

# 读入10铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene 
table(cu_gene %in% rownames(expr.gbm)) 


# 画图
# 热图
library(pheatmap)
expr.cu = expr.gbm[cu_gene,]
n=t(scale(t(expr.cu))) 
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
# 排序
pd.CC = pd.CC[order(pd.CC$cluster),]
ac=data.frame(group=pd.CC$cluster)
rownames(ac)=pd.CC$sample
ac$group
n = n[,rownames(ac)]
pheatmap(n,show_colnames =F,show_rownames = T,
         annotation_col=ac,
         cluster_row = T,cluster_col = F)
pheatmap(n,show_colnames =F,show_rownames = T,
         annotation_col=ac,
         cluster_row = T,cluster_col = F,
         filename = 'Figure3/out_polt/F3_E.pdf')
dev.off()

