# estimate评分

# 设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr) 
library(limma)
library(estimate)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)
library(data.table)



# 读入数据
load("Figure3/out_data/expr.gbm.Rdata")
load("Figure3/out_data/pd.CC.Rdata")
dim(expr.gbm)
expr.gbm[1:4,1:4]
identical(colnames(expr.gbm),pd.CC$sample)

#  do estimate
estimate_RNAseq <- function(RNAseq_logCPM,pro){
  input.f=paste0("Figure4/out_data/",pro,'_estimate_input.txt')
  output.f=paste0("Figure4/out_data/",pro,'_estimate_gene.gct')
  output.ds=paste0("Figure4/out_data/",pro,'_estimate_score.gct')
  write.table(RNAseq_logCPM,file = input.f,sep = '\t',quote = F)
  
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  scores=data.frame(  scores)
  return(scores)
}
pro <- 'CheckCode'
dat <- expr.gbm
scores <- estimate_RNAseq(dat, pro)
head(scores)
save(scores, file = 'Figure4/out_data/estimate_results.Rdata')


# 画图
load('Figure4/out_data/estimate_results.Rdata')
# 肿瘤纯度计算
scores$purity = cos(0.6049872018+0.0001467884 * scores$ESTIMATEScore)
# TIDE计算：
data_tide = read.csv("Input_data/TIDE/export.csv")
rownames(data_tide) = data_tide$Patient
data_tide = data_tide[rownames(scores),]

# # TIS评分
# #  TIS-基因
gene_TIS = c(
  "TIGIT",
  "CD27",
  "CD8A",
  "PDCD1LG2", #(PD-L2)
  "LAG3",
  "CD274", #(PD-L1)
  "CXCR6",
  "CMKLR1",
  "NKG7",
  "CCL5",
  "PSMB10",
  "IDO1",
  "CXCL9",
  "HLA-DQA1",
  "CD276",
  "STAT1",
  "HLA-DRB1",
  "HLA-E")

# 匹配
table(gene_TIS %in% rownames(expr.gbm))
gene_TIS[!gene_TIS %in% rownames(expr.gbm)]
expr_TIS = expr.gbm[gene_TIS,]
# 计算平均值
data_TIS = data.frame(sample = colnames(expr_TIS),
                       TIS_mean = apply(expr_TIS, 2,mean)
) 

# 匹配一下数据
identical(rownames(scores),rownames(data_tide))
data_all = cbind(scores,data_tide)

identical(rownames(data_all),rownames(data_TIS))
data_all = cbind(data_all,data_TIS)
save(data_all,file = "Figure4/out_data/estimate_scores.Rdata")

table(pd.CC$cluster)
identical(row.names(data_all),pd.CC$sample)
data_all$group = pd.CC$cluster

data_all$group=factor(data_all$group, levels=c("C1", "C2"))
group=levels(factor(data_all$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

# 开始画图
colnames(data_all)
score_name = c("StromalScore" , "ImmuneScore" ,  "ESTIMATEScore" ,"purity",
               "TIDE","Dysfunction"  , "Exclusion","TIS_mean")
for(i in score_name){
  # i = "ImmuneScore"
  rt=data_all[,c(i, "group")]
  colnames(rt)=c("IPS", "group")
  gg1=ggviolin(rt, x="group", y="IPS", fill = "group", 
               xlab="group", ylab=i,
               legend.title="group",
               palette=c("#FFFF00", "#FF69B4"),
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons)
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  print(gg1)
  pdf(file=paste0("Figure4/out_polt/TME/","TME_score_",i, ".pdf"), width=6, height=5)
  print(gg1)
  dev.off()
}
