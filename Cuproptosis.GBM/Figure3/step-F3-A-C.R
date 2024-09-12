#目的：CGGA矩阵进行CC分型


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

# 标准化表达矩阵
boxplot(expr.gbm[,1:8],las=2)
max(expr.gbm)
expr.gbm = log2(expr.gbm+1)


# 处理生存数据
rownames(phe.gbm) = phe.gbm$CGGA_ID
phe.gbm = phe.gbm[,c("OS","Censor (alive=0; dead=1)")]
colnames(phe.gbm) = c("time","event")
boxplot(phe.gbm$time)
phe.gbm$time <- phe.gbm$time/365
boxplot(phe.gbm$time)
kp <- phe.gbm$time >0
table(kp)
phe.gbm = na.omit(phe.gbm)
#有几个数据是没有生存数据的，要去除
expr.gbm = expr.gbm[,rownames(phe.gbm)]
dim(expr.gbm)

###########
#保存数据
save(expr.gbm,phe.gbm ,file = 'Figure3/out_data/CGGA_count.Rdata')

######  CC分型

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)

#读入数据
load('Figure3/out_data/CGGA_count.Rdata')
dim(expr.gbm)

# 读入10个铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
expr.cu = expr.gbm[cu_gene,]

# sweep函数减去中位数进行标准化
d = sweep(expr.cu,1, apply(expr.cu,1,median,na.rm=T))
d = as.matrix(expr.cu)

#2、CC分型
dir.create('CC-result')
title="./CC-result"
library(ConsensusClusterPlus)
res <- ConsensusClusterPlus(d, maxK = 7, reps = 1000, pItem = 0.8, 
                            pFeature = 1,corUse = "complete.obs", seed=123456, plot="pdf", writeTable=T,
                            title=title)
# 移动文件
# 手动将文件夹CC-result 移动到 Figure3





