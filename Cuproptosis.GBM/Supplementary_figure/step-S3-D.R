# cu基因的表达

# 设置
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)


# 读入数据
load('Supplementary_figure_3/out_data/count.Rdata')
dim(expr.tumor)
pd.CC = fread( 'Supplementary_figure_3/TCGA-CC-result/TCGA-CC-result.k=2.consensusClass.csv',
               header = F, data.table = F)
pd.CC$cluster = paste0('C',pd.CC$V2)


# 读入10铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
table(cu_gene %in% rownames(expr.tumor)) 

#标准化数据
library(limma)
expr.tumor=normalizeBetweenArrays(expr.tumor)
boxplot(expr.tumor[,1:6],las=2)
exp.cu = expr.tumor[rownames(expr.tumor)  %in% cu_gene,]
data = t(exp.cu)
# 转换长数据
dat <- data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
pd.c1 = pd.CC[pd.CC$cluster == 'C1',]
dat$Group = ifelse(dat$Sample %in% pd.c1$V1,"C1","C2")

library(ggpubr)
library(RColorBrewer)
mypalette <- brewer.pal(3,'Set1')
ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette)+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
ggsave('Supplementary_figure_3/out_polt/S3-D.pdf')