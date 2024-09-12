# H．	密度图展示铁死亡评分high和low两组间免疫检查点的表达差异；

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)


# 读入铁死亡评分high和low两组
# 铜死亡评分
load('Figure6/out_data/cu_score_group.Rdata')
data_cu_score = new_dat

# 读入表达矩阵
load("Figure3/out_data/expr.gbm.Rdata")

# 读入10铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene

# B7-H3 : CD276
# CD279 : PDCD1
# PDCD1Lg2: PDCD1LG2
expr_cu = expr.gbm[cu_gene,]
expr_cu = as.data.frame(t(expr_cu))
new_dat = new_dat[rownames(expr_cu),]
identical(rownames(expr_cu),rownames(new_dat))
expr_cu$cu_group = new_dat$cu_group

# 手动画多个图
i = 10
data = expr_cu[,c(i,11)]
gene_name = colnames(expr_cu)[i]
colnames(data) = c("gene","cu_group")

p10 = ggplot(data, aes(x = gene))+ 
  geom_density(aes(fill = cu_group), alpha=0.4)+
  labs(x = NULL,y= gene_name)+
  theme_bw()+ theme(panel.grid=element_blank())+ 
  xlim(0, 10)+
  scale_fill_manual(values = c( "#FF69B4","#FFFF00"))
p10


p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+plot_layout(ncol = 1,widths = 10,heights = 25)

ggsave("Figure6//out_polt//E_density.pdf")

############################################################################

# 标注p值
#画箱形图
require(tidyr)
require(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tibble)
colour = c("#7CFC00","#DB423E") 

# 宽变长，变成tidy data frame，好用于ggplot2画图
#转换长数据
data_new = melt(expr_cu)
colnames(data_new) = c("group","gene","expression")

#画图
ggplot(data_new,aes(gene,expression,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "gene", y = "expression") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = colour)+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),method = "t.test")


ggsave("Figure6//out_polt//E_density_p-value.pdf")
