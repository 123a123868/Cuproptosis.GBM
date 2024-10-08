# xell计算
# 免疫浸润

devtools::install_github('dviraran/xCell',force = TRUE)

# 设置
library(xCell)

# 读入数据：
# 可以不经过log
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

# 直接计算
xcell_results= xCellAnalysis(expr.gbm)

# 保存数据
save(xcell_results,file = "Figure4/免疫浸润多算法/out_data/xcell_results.Rdata")
