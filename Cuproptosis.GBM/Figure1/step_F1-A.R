# 目的： 比较SCL31A1的差异--cu转运蛋白
# 几组比较：
# TCGA里面的比较
# CGGA的比较
# C1和C2之间的比较
# 高评分和低评分的比较


#设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)
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


# 读入TCGA的GBM数据
library(data.table)
a = fread( 'Input_data/TCGA_GBM/TCGA-GBM.htseq_fpkm.tsv.gz',
           header = T, data.table = F)
dim(a) 
a[1:4,1:4]
#转换一下id  Ensembl_ID 
library(tinyarray)
library(stringi)
expr.TCGA = a
rownames(expr.TCGA) = stri_sub(expr.TCGA$Ensembl_ID,1,15)##保留前15位
expr.TCGA = trans_exp(expr.TCGA)
expr.TCGA = expr.TCGA[,-1]
expr.TCGA[1:4,1:4]

######################
#     SCL31A1
# 表达量的boxp
data = expr.TCGA["SLC31A1",]
data = as.data.frame(t(data)) 
data$group = ifelse(grepl('11A', rownames(data)),'Normal','Tumor') 
colour = c("#7CFC00","#DB423E")
b1 = ggplot(dat = data, aes(group,SLC31A1))+
  geom_boxplot(aes(fill = group))+
  scale_fill_nejm()+
  scale_fill_manual(values= colour)+
  stat_compare_means(method='t.test') + labs(x='expr', title = 'TCGA_SLC31A1')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        legend.position = '')
b1
ggsave(b1,filename = 'Figure1/out_polt/A_1.pdf')

################################################################################################

###  A-2,TCGA_Recurrent_SCL31A1
data_gene = expr.TCGA["SLC31A1",]
data_gene = as.data.frame(t(data_gene)) 
data_gene$group = "na"

# 读入分组数据
pd = read.csv("Input_data/TCGA_GBM//TCGA_N+T_clindata(165).csv")
pd$X = substr(pd$X,1,16)
pd=pd[!duplicated(pd$X),] 
row.names(pd) = pd$X
data_gene = data_gene[row.names(pd),]
identical(row.names(pd),row.names(data_gene))

data_gene$group = pd$Secondary..or.Recurrent
table(data_gene$group)
data_gene$group = ifelse(data_gene$group == "yes","yes" ,"no")
table(data_gene$group)


#### 画图
b2 = ggplot(dat = data_gene, aes(group,SLC31A1))+
  geom_boxplot(aes(fill = group))+
  scale_fill_nejm()+
  scale_fill_manual(values= colour)+
  stat_compare_means(method='t.test') + labs(x='Group', title = 'TCGA_Recurrent_SCL31A1')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        legend.position = '')
b2
ggsave(b2,filename = 'Figure1/out_polt/A_2.pdf')


















