# Chemokines的热图

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)
library(RColorBrewer)
library(tibble)

# 读入数据
load("Figure3/out_data/expr.gbm.Rdata")
load("Figure3/out_data/pd.CC.Rdata")
dim(expr.gbm)
identical(colnames(expr.gbm),pd.CC$sample)

# 画图
# 热图

# 注释
library(pheatmap)
annotation_row = read.csv("Input_data/Chemokines.csv")
row.names(annotation_row) = annotation_row$gene
table(annotation_row$gene %in% rownames(expr.gbm))
rownames(annotation_row)[!rownames(annotation_row) %in% rownames(expr.gbm)]
expr.Chemokines = expr.gbm[rownames(annotation_row),]
expr.Chemokines = na.omit(expr.Chemokines)
rownames(expr.Chemokines)[!rownames(expr.Chemokines) %in% rownames(annotation_row)]
annotation_row = annotation_row[rownames(expr.Chemokines),]
annotation_row = na.omit(annotation_row)
expr.Chemokines = expr.Chemokines[rownames(annotation_row),]
annotation_row = as.data.frame(annotation_row$group)
rownames(annotation_row) = rownames(expr.Chemokines)
colnames(annotation_row) = "group2"
annotation_row$group2

# 画图
pheatmap(expr.Chemokines,show_colnames =F,show_rownames = F) #对那些提取出来的1000个基因所在的每一行取出，组合起来为一个新的表达矩阵
n=t(scale(t(expr.Chemokines))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
#排序一下
pd.CC = pd.CC[order(pd.CC$cluster),]
ac=data.frame(group=pd.CC$cluster)
rownames(ac)=pd.CC$sample
ac$group
n = n[,rownames(ac)]
pheatmap(n,show_colnames =F,show_rownames = T,
         annotation_col=ac,
         annotation_row = annotation_row,
         cluster_row = F,cluster_col = F,
         cellheight = 8.7,
         filename = 'Figure4/out_polt/F4_A.pdf')
dev.off()

# 手动 注释p值
# 箱形图

# 检查数据
dim(expr.Chemokines)
expr.Chemokines[1:4,1:4]
dim(pd.CC)
# 转换为长数据
data = t(expr.Chemokines)
dat <- data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = type,value = Proportion,-Sample)
# 读入分组信息
data_C1 =pd.CC[pd.CC$cluster =="C1",] 
dat$Group = ifelse(dat$Sample %in% data_C1$sample,"C1","C2")
table(dat$Group)

# 画图
mypalette <- brewer.pal(3,'Set1')
ggplot(dat,aes(type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Type", y = "Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = c("#377EB8","#E41A1C"))+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
ggsave('Figure4/out_polt/F4_A_Pvalue.pdf',width = 30)











