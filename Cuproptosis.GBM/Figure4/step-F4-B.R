# HLA家族基因的小提琴图

#设置
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)
library(reshape2)
library(plyr)

# 读入数据
load("Figure3/out_data/expr.gbm.Rdata")
load("Figure3/out_data/pd.CC.Rdata")
dim(expr.gbm)
identical(colnames(expr.gbm),pd.CC$sample)

# 读入HLA家族基因
HLA_gene = c(
  "HLA-B", "HLA-C", "HLA-DOA","HLA-DOB", 
  "HLA-DQB2","HLA-DRB6", "HLA-E", 
  "HLA-F", "HLA-F-AS1", "HLA-G", "HLA-J"
)
HLA_gene
table(HLA_gene %in% rownames(expr.gbm))  
HLA_gene[!HLA_gene %in% rownames(expr.gbm)]
expr.HLA = expr.gbm[HLA_gene,]
expr.HLA = na.omit(expr.HLA)
dim(expr.HLA)
expr.HLA[1:4,1:4]

# 画图
# 小提琴图

# 绘图数据
data = expr.HLA
group=as.data.frame(pd.CC$cluster)
rownames(group) = pd.CC$sample
data_new <- data.frame(t(data))
data_new$sample = row.names(data_new)
data_new <- merge(data_new,group,by.x = "sample",by.y = 0)
# 转换长数据
data_new = melt(data_new)
colnames(data_new) = c("sample","group","gene","expression")
data_new$subject=data_new$sample

# 加载绘图函数
source("Figure4/Function_for_violin_plot.R")

# 绘制小提琴图
library(Rmisc)
Data_summary <- summarySE(
  data_new, measurevar="expression", groupvars=c("group","gene"))
head(Data_summary)
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) 
}
gene_split_violin <- ggplot(data_new,aes(x= gene,y= expression,fill= group))+
  geom_split_violin(trim= F,color="white",scale = "area") + 
  geom_point(data = Data_summary,aes(x= gene, y= expression),pch=19,
             position=position_dodge(0.5),size= 1)+ 
  geom_errorbar(data = Data_summary,aes(ymin = expression-ci, ymax= expression+ci), 
                width= 0.05, 
                position= position_dodge(0.5), 
                color="black",
                alpha = 0.8,
                size= 0.5) +
  scale_fill_manual(values = c("#DB423E","#7CFC00"))+ 
  labs(y=("Log2 expression"),x=NULL,title = "Split violin") + 
  theme_bw()+ mytheme +
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "anova",
                     label.y = max(data_new$expression),
                     hide.ns = T)
gene_split_violin
ggsave(gene_split_violin,
       filename = "Figure4/out_polt/F4_B.pdf",height = 6,width = 8)
