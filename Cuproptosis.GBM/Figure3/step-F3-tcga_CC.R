# TCGA数据集做CC分型


#设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)


# 读入TCGA的GBM数据
library(data.table)
a = fread( 'Input_data/TCGA_GBM/TCGA-GBM.htseq_fpkm.tsv.gz',
           header = T, data.table = F)
dim(a) 
a[1:4,1:4]
#转换一下id  Ensembl_ID 
library(tinyarray)
library(stringi)
expr.TCGA = a
rownames(expr.TCGA) = stri_sub(expr.TCGA$Ensembl_ID,1,15)##保留前15位
expr.TCGA = trans_exp(expr.TCGA)
expr.TCGA = expr.TCGA[,-1]
expr.TCGA[1:4,1:4]

# 临床信息
pd = read.csv("Input_data/TCGA_GBM//TCGA_N+T_clindata(165).csv")
pd$X = substr(pd$X,1,16)
pd=pd[!duplicated(pd$X),] 
row.names(pd) = pd$X
colnames(pd)
pd = pd[c("Secondary..or.Recurrent",
           "Age.at..Procedure",
           "Vital.Status",
           "OS")]
colnames(pd) = c("Recurrent","Age","event","os.time")
pd$os.time
boxplot(pd$os.time)
pd = dplyr::filter(pd, !is.na(os.time))
pd$os.time = pd$os.time/365
boxplot(pd$os.time)
save(pd,file = "Figure3/out_data/TCGA_GBM_pd.Rdata")


# 匹配一下
expr.TCGA = expr.TCGA[,rownames(pd)]
identical(rownames(pd),colnames(expr.TCGA))

# 对TCGA数据集做CC分型
# 标准化表达矩阵
boxplot(expr.TCGA[,1:8],las=2)
max(expr.TCGA)

# 读入10个铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
expr.cu = expr.TCGA[cu_gene,]

# sweep函数减去中位数进行标准化
d = sweep(expr.cu,1, apply(expr.cu,1,median,na.rm=T))
d = as.matrix(expr.cu)

#2、CC分型
dir.create('TCGA_gbm_CC_result')
title="./TCGA_gbm_CC_result"
library(ConsensusClusterPlus)
res <- ConsensusClusterPlus(d, maxK = 7, reps = 1000, pItem = 0.8, 
                            pFeature = 1,corUse = "complete.obs", seed=123456, plot="pdf", writeTable=T,
                            title=title)
# 移动文件
# 手动将文件夹CC-result 移动到 Figure3



