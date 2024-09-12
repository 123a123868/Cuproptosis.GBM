# 铜死亡CGGA(训练集)和TCGA(验证集)麻烦你再补一个这个森林图，
# 因素就选cluster(C1,C2)、Age(>=60, <60)、PRS type(原发，复发)。


#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
library(finalfit) # 关键包
library(ggplot2)

# 先画铜死亡CGGA(训练集)
#读入CC分型信息
pd.CC = fread( 'Figure3/CC-result/CC-result.k=2.consensusClass.csv',
               header = F, data.table = F)
colnames(pd.CC) = c('sample','cluster')
pd.CC$cluster = paste0('C',pd.CC$cluster) 
rownames(pd.CC) = pd.CC$sample

# 临床信息
phe = fread("Input_data/CGGA_GBM/CGGA.mRNAseq_693_clinical.20200506.txt/CGGA.mRNAseq_693_clinical.20200506.txt",
              header = T, sep = '\t',data.table = F)
phe[1:4,1:4]
rownames(phe) = phe$CGGA_ID

# 匹配
phe = phe[pd.CC$sample,]
identical(phe$CGGA_ID,pd.CC$sample)
phe$cc_cluster = pd.CC$cluster
colnames(phe)
pd = phe[c("PRS_type",
          "Gender",
          "OS",
          "Age",
          "Censor (alive=0; dead=1)",
          "cc_cluster")]

#生存信息
boxplot(pd$OS)  #年的时间
pd$time <- pd$OS/365
head(pd)
boxplot(pd$time) 
colnames(pd)[5] = "event"
pd$Age = ifelse(pd$Age <65,"young","old")


# 开始画图
colnames(pd)
explanatory = c("PRS_type","Age","cc_cluster")
dependent = 'event'
data_forest = pd
data_forest %>%
  or_plot(dependent, explanatory)
pdf("Figure4/out_polt/forest_CGGA.pdf",width = 10,height=5)
data_forest %>%
  or_plot(dependent, explanatory)
dev.off()



########################################################################################\
#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
library(finalfit) # 关键包
library(ggplot2)

# TCGA(验证集)
#读入CC分型信息
pd.CC = fread( 'Figure3/TCGA_gbm_CC_result/TCGA_gbm_CC_result.k=2.consensusClass.csv',
               header = F, data.table = F)
colnames(pd.CC) = c('sample','cluster')
pd.CC$cluster = paste0('C',pd.CC$cluster) 
rownames(pd.CC) = pd.CC$sample

# 读入临床信息
load("Figure3/out_data/TCGA_GBM_pd.Rdata")
pd$Age = ifelse(pd$Age <65,"young","old")
pd$event = ifelse(pd$event == "DECEASED","1","0")
class(pd$event)
pd$event = as.numeric(pd$event)
identical(rownames(pd),rownames(pd.CC))
pd$cc_cluster = pd.CC$cluster
pd = na.omit(pd)
table(pd$Recurrent)
pd$Recurrent = ifelse(pd$Recurrent == "yes","yes","No")

# 开始画图
colnames(pd)
explanatory = c("Recurrent","Age","cc_cluster")
dependent = 'event'
data_forest = pd
data_forest %>%
  or_plot(dependent, explanatory)
pdf("Figure4/out_polt/forest_TCGA.pdf",width = 10,height=5)
data_forest %>%
  or_plot(dependent, explanatory)
dev.off()








