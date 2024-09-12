# 免疫检查点的boxpolt展示

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
library(survminer) 
library(survival)
library(loose.rock)
library(futile.logger) 
library(glmSparseNet)
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)

# 读入表达矩阵
load("Figure3/out_data/expr.gbm.Rdata")
dim(expr.gbm)

# 读入铁死亡评分high和low两组
load("Figure6/out_data/cu_score_group.Rdata")
head(new_dat)

# 免疫检查点基因
immune_gene = c("CTLA4","CD276","PDCD1","CD274","HAVCR2","CD27",
                "CD40","CD47","CD70","ICOS","IDO1","LAG3","LGALS9",
                "PDCD1","PDCD1LG2","PVR","VTCN1","TNFSF14","TNFSF18")

table(immune_gene %in% rownames(expr.gbm))  #都能在里面找到表达。
immune_gene[!immune_gene %in% rownames(expr.gbm)]
expr.immune = expr.gbm[immune_gene,]


#画箱形图
require(tidyr)
require(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tibble)
mypalette <- brewer.pal(3,'Set1')

# 宽变长，变成tidy data frame，好用于ggplot2画图
expr.immune = as.data.frame(t(expr.immune))
new_dat = new_dat[rownames(expr.immune),]
identical(rownames(expr.immune),rownames(new_dat))
expr.immune$cu_group = new_dat$cu_group

#转换长数据
data_new = melt(expr.immune)
colnames(data_new) = c("group","gene","expression")


#画图
mypalette <-c("#FF69B4","#FFFF00")
ggplot(data_new,aes(gene,expression,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "gene", y = "expression") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = c("#FF69B4","#FFFF00"))+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),method = "t.test")
ggsave('Figure6/out_polt//boxplot_immune_point.pdf')


