# tcga的数据进行CC分型
#设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)

# 读入数据
library(data.table)
a = fread( 'Input_data/TCGA_GBM/TCGA-GBM.htseq_counts.tsv.gz',
           header = T, data.table = F)
dim(a) 
a[1:4,1:4]
# 转换id  
library(tinyarray)
library(stringi)
expr = a
rownames(expr) = stri_sub(expr$Ensembl_ID,1,15)
expr = trans_exp(expr) 
expr = expr[,-1]
expr[1:4,1:4]
## tumor和normal的分组 
dim(expr)
colnames(expr)
group = ifelse(grepl('11A', colnames(expr)),'Normal','Tumor') 
table(group)
expr.tumor = expr[,group == 'Tumor']

# 读入10铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
table(cu_gene %in% rownames(expr.tumor)) 
expr.cu = expr.tumor[cu_gene,]
save(expr.tumor,cu_gene,file = 'Supplementary_figure_3//out_data/count.Rdata')
# sweep函数减去中位数进行标准化
d = sweep(expr.cu,1, apply(expr.cu,1,median,na.rm=T))
d = as.matrix(expr.cu)
# CC分型
setwd("Supplementary_figure_3/")
dir.create('TCGA-CC-result')
title="./TCGA-CC-result"
library(ConsensusClusterPlus)
res <- ConsensusClusterPlus(d, maxK = 7, reps = 1000, pItem = 0.8, 
                            pFeature = 1,corUse = "complete.obs", seed=123456, plot="pdf", writeTable=T,
                            title=title)
setwd("../")
