## 免疫检查点的箱形图

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)

# 读入数据
load("Figure3/out_data/expr.gbm.Rdata")
load("Figure3/out_data/pd.CC.Rdata")
dim(expr.gbm)
identical(colnames(expr.gbm),pd.CC$sample)

# 读入免疫检查点基因
ICG_gene = c("CTLA4", "CD27", "LGALS9", "LAG3", "CD274", "PDCD1", "VTCN1")
table(ICG_gene %in% rownames(expr.gbm))  #都能在里面找到表达。
ICG_gene[!ICG_gene %in% rownames(expr.gbm)]
expr.ICG = expr.gbm[ICG_gene,]

# 绘图
# 箱形图

# 绘图数据
# 转化长数据
data = t(expr.ICG)
dat <- data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = type,value = Proportion,-Sample)
#读入分组信息
data_C1 =pd.CC[pd.CC$cluster =="C1",] 
dat$Group = ifelse(dat$Sample %in% data_C1$sample,"C1","C2")
table(dat$Group)

# 画图
ggplot(dat,aes(type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Type", y = "Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = c("#DB423E","#7CFC00"))+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
ggsave('Figure4/out_polt/F4_E.pdf')




