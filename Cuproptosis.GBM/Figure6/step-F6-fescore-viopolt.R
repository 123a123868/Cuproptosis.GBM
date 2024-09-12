# 图C：C．	FUSCC C1和C2铁死亡分数比较（小提琴图）；

# 设置
rm(list = ls())  
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

# 铁死亡评分
load('Figure6/out_data/cu_score_group.Rdata')

# 读入CC分型信息
pd.CC = fread( 'Figure3/CC-result/CC-result.k=2.consensusClass.csv',
               header = F, data.table = F)
colnames(pd.CC) = c('sample','cluster')
pd.CC$cluster = paste0('C',pd.CC$cluster) 

# 匹配一下
rownames(pd.CC) = pd.CC$sample
new_dat = new_dat[pd.CC$sample,]
identical(rownames(pd.CC),rownames(new_dat))
pd.CC$cu_score = new_dat$cu_score

# 画小提琴图
data.score  = pd.CC
data.score$group = data.score$cluster
data.score$group=factor(data.score$group, levels=c("C1", "C2"))
group=levels(factor(data.score$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
colnames(data.score)

# 开始画图
# c("#FFFF00", "#FF69B4")，即C1用黄色，C2用粉色哈，不然好多AI没法换色
for(i in colnames(data.score)[3]){
  # i = "cu_score"
  rt=data.score[,c(i, "group")]
  colnames(rt)=c("IPS", "group")
  gg1=ggviolin(rt, x="group", y="IPS", fill = "group", 
               xlab="group", ylab=i,
               legend.title="group",
               palette=c("#FFFF00", "#FF69B4"),
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons)
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  print(gg1)
  
  pdf(file=paste0("Figure6/out_polt/","C_cu_score_CC", ".pdf"), width=6, height=5)
  print(gg1)
  dev.off()
}
