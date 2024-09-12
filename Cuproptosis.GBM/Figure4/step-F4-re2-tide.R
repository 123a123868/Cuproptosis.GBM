# TIDE计算

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)


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

# 临床信息
phe.a = fread("Input_data/CGGA_GBM/CGGA.mRNAseq_693_clinical.20200506.txt/CGGA.mRNAseq_693_clinical.20200506.txt",
              header = T, sep = '\t',data.table = F)
phe.a[1:4,1:4]
phe.a = phe.a[!phe.a$CGGA_ID %in% 'CGGA_1615',]  #该样本特殊

#筛选出GBM数据
library(dplyr)
phe.gbm = filter(phe.a,Histology== 'GBM' |Histology == 'rGBM')
table(phe.gbm$Histology)
expr.gbm = expr.a[,phe.gbm$CGGA_ID]
dim(expr.gbm)
expr.gbm[1:4,1:4]
dim(expr.gbm)


# 标准化数据
library(preprocessCore)
expr.gbm = as.matrix(expr.gbm)
expr.gbm = normalize.quantiles(expr.gbm)
boxplot(expr.gbm[,1:8],las=2)
max(expr.gbm)
min(expr.gbm)

# 所有基因的表达值都是通过减去数据集中所有样本的平均值来规范化的
# 直接减去平均值
library(preprocessCore)
data = as.matrix(expr.gbm)
data = normalize.quantiles(data)
boxplot(data[,1:8],las=2)
colnames(data) = colnames(expr.gbm)
rownames(data) = rownames(expr.gbm)
data[1:4,1:4]
mean(data)
data = data-mean(data)
data[1:4,1:4]
boxplot(data[,1:8],las=2)
max(data)
min(data)
mean(data)
write.table(data,sep ="\t",file = "Figure4/out_data/expr.gbm_scale.txt",quote=F)









