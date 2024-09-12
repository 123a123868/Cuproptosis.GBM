#TCGA队列 免疫检查位点


#读入表达矩阵
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)

# 读入数据
load('Supplementary_figure_3/out_data/count.Rdata')
dim(expr.tumor)
pd.CC = fread( 'Supplementary_figure_3/TCGA-CC-result/TCGA-CC-result.k=2.consensusClass.csv',
               header = F, data.table = F)
colnames(pd.CC) = c('sample','cluster')
pd.CC$cluster = paste0('C',pd.CC$cluster) 
identical(pd.CC$sample,colnames(expr.tumor))
gp = pd.CC$cluster
table(gp)

#读入ICD基因
gene.ICB = c('CD47','CD70','IDO1','LAG3','PDCD1','VTCN1','CD27', 'LGALS9')
gene.ICB = c("CTLA4", 
             "CD276", 
             "CD274", 
             "HAVCR2", 
             "CD27", 
             "CD40", 
             "CD47", 
             "CD70", 
             "ICOS", 
             "IDO1", 
             "LAG3", 
             "LGALS9", 
             "PDCD1", 
             "PDCD1LG2", 
             "PVR", 
             "VTCN1", 
             "TNFSF14", 
             "TNFSF18")
table(gene.ICB %in% rownames(expr.tumor))
gene.ICB[!gene.ICB %in% rownames(expr.tumor)]

#取基因的表达矩阵
exp.ICB = expr.tumor[rownames(expr.tumor)  %in% gene.ICB,]
# rownames(exp.ICB) = c('PD_L1','CTLA4','IOD1','TIM3',"LAG3","Gal9" ,'PD1','PD-L2')
# exp.ICB = exp.ICB[!rownames(exp.ICB) %in% c("IOD1","PD-L2") ,]


data = t(exp.ICB)
rownames(exp.ICB)

# 画图
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)


dat <- data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

#读入CC分型信息
pd.c1 = pd.CC[pd.CC$cluster == 'C1',]
dat$Group = ifelse(dat$Sample %in% pd.c1$sample,"C1","C2")


library(ggpubr)
library(RColorBrewer)
mypalette <- brewer.pal(3,'Set1')
ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Expression") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette)+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "wilcox.test")
ggsave('Supplementary_figure_3/out_polt//tcga_CC_PD1.pdf')







