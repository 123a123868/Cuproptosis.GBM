# IMvigor210CoreBiologies

#设置
rm(list = ls()) 
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
#install.packages("12-IMvigor210CoreBiologies/IMvigor210CoreBiologies_1.0.0.tar.gz", repos=NULL)
# if (!requireNamespace("IOBR", quietly = TRUE))
#   devtools::install_github("IOBR/IOBR")
library(IOBR) 
library(IMvigor210CoreBiologies)
data(cds)
expMatrix <- counts(cds)
eff_length2 <- fData(cds)[,c("entrez_id","length","symbol")]
rownames(eff_length2) <- eff_length2$entrez_id
head(eff_length2)
feature_ids <- rownames(expMatrix)
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2),]
mm <- match(rownames(expMatrix),rownames(eff_length2))
eff_length2 <- eff_length2[mm,]

x <- expMatrix/eff_length2$length
eset <- t(t(x)/colSums(x))*1e6
summary(duplicated(rownames(eset)))

eset <- IOBR::anno_eset(eset = eset,
                        annotation = eff_length2,
                        symbol = "symbol",
                        probe = "entrez_id",
                        method = "mean")
tumor_type <- "blca"
if(max(eset)>100) eset <- log2(eset+1)
dim(eset)


######
#  表型文件
pdata <- pData(cds)
colnames(pdata)
table(pdata$`Best Confirmed Overall Response`)
colnames(pdata) <- gsub(colnames(pdata),pattern = " ",replacement = "_")
colnames(pdata) 
# 修改列名：
names(pdata)[names(pdata) == 'censOS'] <- 'event'
names(pdata)[names(pdata) == 'os'] <- 'time'
names(pdata)[names(pdata) == 'binaryResponse'] <- 'BOR_binary'
pdata<-pdata[!is.na(pdata$BOR_binary),]
pdata$pd1<-ifelse(pdata$BOR_binary=="CR/PR","Responder","NonResponder")
pdata$ID = rownames(pdata)
table(pdata$pd1)
table(pdata$BOR_binary)
table(colnames(eset) %in% pdata$ID)
expr = eset[,pdata$ID]
dim(expr)
save(expr,pdata,file = "Figure7/out_data//IMvigor_count.Rdata")

########################################################################
# cu评分
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


load("Figure7/out_data//IMvigor_count.Rdata")
dim(expr)
expr[1:4,1:4]
dim(pdata)
pdata[1:4,1:4]
rownames(pdata) = pdata$ID
identical(colnames(expr), rownames(pdata))

colnames(pdata)

# cu评分
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
kp <- cu_gene
cdat <- expr[kp, ]
xdata <- t(cdat)
##    计算cu评分。用e，coefficient的评分公式
xdata = as.data.frame(xdata)
xdata$cu_score = exp(-0.262490103*xdata$CDKN2A
                          -0.266355127*xdata$DLAT
                          +0.184258117*xdata$DLD
                          -0.374608653*xdata$FDX1
                          -0.086511329*xdata$GLS
                          -0.242718907*xdata$LIAS
                          +0.508961778*xdata$LIPT1
                          +0.220024283*xdata$MTF1
                          +0.001332578*xdata$PDHA1
                          +0.149181298*xdata$PDHB)
ydata <- pdata[,c("time","event")]
head(ydata)
identical(rownames(xdata), rownames(ydata))
## time transform to year
head(ydata)
boxplot(ydata$time)
ydata$time = ydata$time/12
kp <- ydata$time >0
table(kp)
xdata <- xdata[kp,]
ydata <- ydata[kp,]

#   合并xdata和ydata
identical(rownames(xdata),rownames(ydata))
dat_cox <- cbind(ydata,xdata)
head(dat_cox)

##  surv_cutpoint函数，高低分组的最佳切点
res.cut <- surv_cutpoint(dat_cox, 
                         time = "time", 
                         event = "event", 
                         variables = "cu_score" 
)
summary(res.cut) #查看数据最佳截断点及统计量
cu_group <- as.factor(ifelse(dat_cox$cu_score > 0.6946737,'High','Low'))
dat_cox$cu_group <- cu_group 
table(dat_cox$cu_group)
save(dat_cox,file = "Figure7/out_data//cu_scroe_group.Rdata")






























