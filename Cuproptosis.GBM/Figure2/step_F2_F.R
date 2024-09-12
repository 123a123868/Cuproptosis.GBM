# 目的
# 10个铜基因在TGCA-GBM的表达

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)
library(ggplot2)
library(ggsci)
library(ggstatsplot)

#####################
#1、处理数据
#    TCGA
#####################
#读入数据
library(data.table)
a = fread( 'Input_data/TCGA_GBM/TCGA-GBM.htseq_fpkm.tsv.gz',
           header = T, data.table = F)
dim(a) 
a[1:4,1:4]
# 转换一下id  
library(tinyarray)
library(stringi)
expr.TCGA = a
rownames(expr.TCGA) = stri_sub(expr.TCGA$Ensembl_ID,1,15) 
expr.TCGA = trans_exp(expr.TCGA)
expr.TCGA = expr.TCGA[,-1]
expr.TCGA[1:4,1:4]

# 读入10铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
table(cu_gene %in% rownames(expr.TCGA))  
boxplot(expr.TCGA[,1:8],las=2)
expr.TCGA = expr.TCGA[cu_gene,]


#画箱形图
require(tidyr)
require(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tibble)
colour = c("#7CFC00","#DB423E") 

# 宽变长，变成tidy data frame，好用于ggplot2画图
data = t(expr.TCGA)
dat <- data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = type,value = Proportion,-Sample)
#读入分组信息
dat$Group = ifelse(grepl('11A', dat$Sample),'Normal','Tumor') 
#画图
ggplot(dat,aes(type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Type", y = "Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = colour)+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "t.test")
ggsave('Figure2/out_polt/F2_F.pdf')
