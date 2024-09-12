#制作输入文件
# 需要的是nolog的数据
# 读入数据
#设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)
library(tinyarray)
library(tidyverse)

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
expr_ciber = rownames_to_column(expr.gbm)


###### 开始计算
write.table(expr_ciber,file = "Figure4/免疫浸润多算法/input_data/CIBERSORT//exp.txt",row.names = F,quote = F,sep = "\t")

source("Figure4/免疫浸润多算法/input_data/CIBERSORT//CIBERSORT_2.R")

if(F){
  TME.results = CIBERSORT("Figure4/免疫浸润多算法/input_data/CIBERSORT//LM22-ref.txt", 
                          "Figure4/免疫浸润多算法/input_data/CIBERSORT//exp.txt" , 
                          perm = 1000, 
                          QN = T)
  ciber.results = TME.results
  save(ciber.results,file = "Figure4/免疫浸润多算法/out_data/ciber.results.Rdata")
}
load("Figure4/免疫浸润多算法/input_data/CIBERSORT//ciber_CHOL.Rdata")
TME.results[1:4,1:4]
re <- TME.results[,-(23:25)]
