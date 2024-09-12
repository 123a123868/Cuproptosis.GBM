# 药物敏感性分析

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

# 读入GDSC数据
dir= "Input_data/oncoPredict_GDSC/Training Data/"
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

# 表达矩阵和分组
load("Figure3/out_data/expr.gbm.Rdata")
load("Figure3/out_data/pd.CC.Rdata")
dim(expr.gbm)
identical(colnames(expr.gbm),pd.CC$sample)

# 药敏分析
testExpr = as.matrix(expr.gbm)
setwd("Figure5/out_data/")
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData')
setwd("../../")
getwd()

#批量计算差异
library(foreign)
library(car) 
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

# 读入药敏数据
data_drug = read.csv("Figure5/out_data/calcPhenotype_Output/DrugPredictions.csv")
dim(data_drug)
data_drug[1:4,1:4]
rownames(data_drug) = data_drug$X
data_drug = data_drug[,-1]
data_drug[1:4,1:4]
data_drug = data_drug[pd.CC$sample,]
identical(pd.CC$sample,rownames(data_drug))
data_drug$group = pd.CC$cluster

# 计算p值，
dim(data_drug)
pval_pearson <- apply(data_drug[1:198], 2, function(x) 
  (wilcox.test(formula = x ~ group, data = data_drug))[3])
library (plyr)
df <- ldply (pval_pearson, data.frame)

options(digits = 4)
df$p.value
df$p= round(df$p.value, 3)
table(df$p < 0.05)
write.csv(df,file = 'Figure5/out_data/drug_CC_p.value.csv')





