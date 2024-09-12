

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)


# 10个铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene


