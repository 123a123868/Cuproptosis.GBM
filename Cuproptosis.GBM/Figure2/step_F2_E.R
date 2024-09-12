# 目的
# 10个铜基因在CGGA的表达
# 合并GTEX 数据集和CGGA的数据集


#设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)


#####################
#1、处理数据
#    CGGA数据库
#####################
#读入数据
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

#筛选出GBM数据
library(dplyr)
phe.gbm = filter(phe.a,Histology== 'GBM' |Histology == 'rGBM')
table(phe.gbm$Histology)
expr.gbm = expr.a[,phe.gbm$CGGA_ID]
dim(expr.gbm)
expr.gbm[1:4,1:4]
dim(expr.gbm)

#####################
#2、处理数据
#   GTEX数据库
#####################
#读入GTEx数据
expr.gtex=fread("Input_data/GTEX//gtex_RSEM_gene_fpkm.gz",header = T, sep = '\t',data.table = F)
dim(expr.gtex)
expr.gtex[1:4,1:4]

#转换一下id  Ensembl_ID---
library(tinyarray)
library(stringi)
rownames(expr.gtex) = stri_sub(expr.gtex$sample,1,15)##保留前15位
expr.gtex = expr.gtex[,-1]
expr.gtex[1:3,1:3]


#读取gtex的临床样本注释信息
#选取Thyroid的样本。366个
phe.gtex = fread("Input_data/GTEX//GTEX_phenotype.gz",header = T, sep = '\t',data.table = F)
colnames(phe.gtex)
colnames(phe.gtex)=c("Sample","body_site","primary_site","gender","patient","cohort")
colnames(phe.gtex)
rownames(phe.gtex)=phe.gtex$Sample
table(phe.gtex$primary_site)
table(phe.gtex$body_site)
phe.gtex.Brain <- phe.gtex[c(grep('*Brain - Cortex',phe.gtex$body_site),
                             grep('*Brain - Frontal Cortex',phe.gtex$body_site)),]
table(phe.gtex.Brain$body_site)


#筛选Thyroid的样本的表达矩阵
table(colnames(expr.gtex) %in% phe.gtex.Brain$Sample)
expr.gtex = expr.gtex[,colnames(expr.gtex) %in% phe.gtex.Brain$Sample]
expr.gtex[1:4,1:4]
phe.gtex.Brain2 = phe.gtex.Brain[colnames(expr.gtex),]


#    两个数据集合并在一起
#    我们从官网可以看到gtex是按照log2(fpkm+0.001)处理的  
#    CGGA的是FPKM 的，且完全没有经过处理的
#   处理一下表达矩阵。标准化
boxplot(expr.gbm[,1:8],las=2)
boxplot(expr.gtex[,1:8],las=2)

#①、逆log
expr.gtex2 = 2**expr.gtex-0.001  
boxplot(expr.gtex2[,1:8],las=2)  #变为原始的fpkm。

# 两个矩阵合并起来。
expr.gbm[1:4,1:4]
expr.gtex2[1:4,1:4]

library(clusterProfiler)
library(org.Hs.eg.db)
genenames<-bitr(rownames(expr.gtex2),fromType = "ENSEMBL",
                toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
genenames=genenames[!duplicated(genenames$SYMBOL),] 
genenames=genenames[!duplicated(genenames$ENSEMBL),] 
expr.gtex2 = expr.gtex2[genenames$ENSEMBL,]
identical(rownames(expr.gtex2),genenames$ENSEMBL)
rownames(expr.gtex2) = genenames$SYMBOL

kp = intersect(rownames(expr.gbm),rownames(expr.gtex2))
expr.gbm = expr.gbm[kp,]
expr.gtex2 = expr.gtex2[kp,]
identical(rownames(expr.gbm),rownames(expr.gtex2))

expr.all = cbind(expr.gbm,expr.gtex2)
dim(expr.all)
dim(expr.gbm)
dim(expr.gtex2)
expr.all[1:4,1:4]

# 统一进行log+1
expr.all = log2(expr.all+1)
boxplot(expr.all[,1:8],las=2) 
boxplot(expr.all[,450:456],las=2) 


# 分组信息
colnames(expr.all)
group = ifelse(grepl('CGGA', colnames(expr.all)),'Tumor','Normal') 
table(group)


######################################################
#铜基因在TCGA里面的表达情况
# 10个铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
table(cu_gene %in% rownames(expr.all))
boxplot(expr.all[,1:8],las=2)
expr.all = expr.all[cu_gene,]


#画箱形图
require(tidyr)
require(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tibble)
colour = c("#7CFC00","#DB423E") 

# 宽变长，变成tidy data frame，好用于ggplot2画图
data = t(expr.all)
dat <- data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = type,value = Proportion,-Sample)
#读入分组信息
dat$Group = ifelse(grepl('CGGA', dat$Sample),'Tumor','Normal') 
#画图
ggplot(dat,aes(type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Type", y = "Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = colour)+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "t.test")
ggsave('Figure2/out_polt/F2_E.pdf')
