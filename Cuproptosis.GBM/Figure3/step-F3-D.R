# C1和C2的生存分析

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
library(survminer) 
library(survival)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)

# 读入CC分型信息
pd.CC = fread( 'Figure3/CC-result/CC-result.k=2.consensusClass.csv',
               header = F, data.table = F)
colnames(pd.CC) = c('sample','cluster')
pd.CC$cluster = paste0('C',pd.CC$cluster) 

# 读入临床信息
phe.a = fread("Input_data/CGGA_GBM/CGGA.mRNAseq_693_clinical.20200506.txt/CGGA.mRNAseq_693_clinical.20200506.txt",
              header = T, sep = '\t',data.table = F)
phe.a[1:4,1:4]
rownames(phe.a) = phe.a$CGGA_ID
phe.a = phe.a[,c("OS","Censor (alive=0; dead=1)")]
colnames(phe.a) = c("time","event")
boxplot(phe.a$time)
phe.a$time <- phe.a$time/365
boxplot(phe.a$time)
kp <- phe.a$time >0
table(kp)
phe.a = na.omit(phe.a)

# 匹配取交集
phe.gbm = phe.a[pd.CC$sample,]
identical(rownames(phe.gbm),pd.CC$sample)
phe.gbm$cluster = pd.CC$cluster
phe.gbm[1:4,1:3]

# 画图
#生存曲线
colnames(phe.gbm)
sfit <- survfit(Surv(time, event)~cluster, data=phe.gbm)
summary(sfit)
survp=ggsurvplot(
  sfit,                 
  legend.title = 'cluster', 
  pval = T, 
  risk.table = TRUE, 
  risk.table.y.text = F,
  xlab = "Time in years", 
  break.time.by = 1,
  size = 1.5, 
  ggtheme = theme_ggstatsplot(),
  palette="nejm", 
)
print(survp)
pdf('Figure3/out_polt/F3_D.pdf' , onefile = F)
print(survp)
dev.off()













