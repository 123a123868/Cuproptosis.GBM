# 建模

# 设置
rm(list=ls())
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


# 读入数据
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
#处理一下表达矩阵。标准化
boxplot(expr.gbm[,1:8],las=2)
expr.gbm = log2(expr.gbm+1)

# 筛选生存时间不是0
rownames(phe.gbm) = phe.gbm$CGGA_ID
phe.gbm = phe.gbm[,c("OS","Censor (alive=0; dead=1)")]
colnames(phe.gbm) = c("time","event")
boxplot(phe.gbm$time)
phe.gbm$time <- phe.gbm$time/365
boxplot(phe.gbm$time)
kp <- phe.gbm$time >0
table(kp)
phe.gbm = na.omit(phe.gbm)
# 有几个样本是没有生存数据
expr.gbm = expr.gbm[,rownames(phe.gbm)]
dim(expr.gbm)
save(expr.gbm,phe.gbm ,file = 'Figure6/out_data/lasso_count.Rdata')


# batch km by cu_gene

# 读入数据
rm(list = ls()) 
# 读入10铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
# 读入矩阵和生存信息
load(file = 'Figure6/out_data/lasso_count.Rdata')
expr.gbm[1:4,1:4] 
head(phe.gbm)
boxplot(phe.gbm$time)
identical(colnames(expr.gbm), rownames(phe.gbm))

# 2. lapply km plot----
model_list <- lapply(cu_gene, function(i){
  # i = cu_gene[1]
  survival_dat = phe.gbm
  gene = as.numeric(expr.gbm[i,])
  survival_dat$gene = ifelse(gene > median(gene),'high','low')
  table(survival_dat$gene)
  library(survival)
  fit <- survfit(Surv(time, event) ~ gene,
                 data = survival_dat)
  survp = ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
                     legend.title = i, #定义图例的名称 
                     # legend = "top",#图例位置
                     # legend.labs = c('High', 'Low'),
                     pval = T, #在图上添加log rank检验的p值
                     # pval.method = TRUE,#添加p值的检验方法
                     risk.table = TRUE, 
                     risk.table.y.text = F,
                     xlab = "Time in years", #x轴标题
                     # xlim = c(0, 10), #展示x轴的范围
                     break.time.by = 1, #x轴间隔
                     size = 1.5, #线条大小
                     ggtheme = theme_ggstatsplot(),
                     palette="nejm", #配色
  )
  return(survp)               
}) 

x=3;y=4
all_plot <- arrange_ggsurvplots(model_list,print = F,ncol =x, nrow = y,
                                risk.table.height = 0.3,
                                surv.plot.height = 0.7)
# all_plot  
x=30;y=30
ggsave(all_plot,filename = 'Figure6/out_polt/10cu_genes_KM.pdf' ,
       height = x, width = y )

# 选择3个生存差异的基因
cu_gene = c("CDKN2A","LIAS","DLD")

#####################################################################
# lasso

# 读入数据
rm(list = ls()) 
# 选择3个生存差异的基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene

cu_gene = c("CDKN2A","LIAS","DLD")
# 读入矩阵和生存信息
load(file = 'Figure6/out_data/lasso_count.Rdata')
expr.gbm[1:4,1:4] 
head(phe.gbm)
boxplot(phe.gbm$time)
identical(colnames(expr.gbm), rownames(phe.gbm))

# tidy data for LASSO----
kp <- cu_gene
cdat <- expr.gbm[kp, ]
xdata <- t(cdat)
ydata <- phe.gbm
head(ydata)
identical(rownames(xdata), rownames(ydata))
## time transform to year
head(ydata)
boxplot(ydata$time)
kp <- ydata$time >0
table(kp)
xdata <- xdata[kp,]
ydata <- ydata[kp,]

# run lasso model
source('Figure6//run_LASSO_cox.R') 
lasso_cox(xdata,
          ydata)
save(xdata,ydata,file = 'Figure6/out_data//input_model_data.Rdata')
gene_model <- read.table('Figure6/out_data//model_lasso_cg_genes.txt')[, 1]
gene_model   # lasso筛选 后的基因列表




