# cu评分分组的细胞浸润

# 设置
rm(list = ls())
options(stringsAsFactors = F)
library(pheatmap) 
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)
library(data.table)
library(GSVA)
library(limma)
library(GSEABase)

# 读入数据
load('Figure6/out_data/lasso_count.Rdata')
dim(expr.gbm)
expr.gbm[1:4,1:4]
# 读入免疫细胞浸润数据
gmtFile <- "Input_data/Immune_cell_infiltration//immune.gmt"   
dimnames <- list(rownames(expr.gbm),colnames(expr.gbm))
mat <- matrix(as.numeric(as.matrix(expr.gbm)),nrow=nrow(expr.gbm),dimnames=dimnames)
mat <- avereps(mat)
mat <- mat[rowMeans(mat)>0,]
geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())

# do gsva
ssgseaScore <- gsva(mat, geneSet, 
                    method='ssgsea', 
                    kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaOut <- normalize(ssgseaScore)
ssgseaOut <- rbind(id=colnames(ssgseaOut),ssgseaOut)
dim(ssgseaOut)
rownames(ssgseaOut)
ssgseaOut[1:4,1:4]

# save data----
write.table(ssgseaOut,file="Figure7/out_data/ssgseaOut.txt",
            sep="\t",quote=F,col.names=F)


################################################################################
# 可视化画图
library(pheatmap) 
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)
library(data.table)

#目的：高低风险的免疫浸润差异。
#设置
#目的：可视化gsva的结果，出图
#数据来源：
#保存数据：

rm(list = ls())  
# 0. 倒入数据
input <- "Figure7/out_data/ssgseaOut.txt"  
immune <- read.table( input ,sep="\t",header=T,row.names=1,check.names=F)
immune[1:4, 1:4]

#读入高低分型信息
load("Figure6/out_data/cu_score_group.Rdata")
pd.cu = new_dat
pd.cu$sample = rownames(pd.cu)
pd.cu$cluster  = pd.cu$cu_group

# 1. overlap----
pd.cu = pd.cu[order(pd.cu$cluster),]
immune = immune[,pd.cu$sample]
Type <- pd.cu$cluster
names(Type) <- colnames(immune)
Type <- as.data.frame(Type)

# 3. vioplot----
outTab=data.frame()
table(Type$Type)
High = 87
Low = 150
immune <- data.frame(t(immune),check.names = F)

pdf("Figure7/out_polt/F7_A.pdf",height=7,width=15)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(immune))
y=c(1:ncol(immune))
plot(x,y,
     xlim=c(0,72),ylim=c(min(immune),max(immune)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=24,
     col="white",
     xaxt="n")

library(vioplot)
pFilter=0.05
for(i in 1:ncol(immune)){
  # i = 1
  if(sd(immune[1:High,i])==0){
    immune[1,i]=0.001
  }
  if(sd(immune[(High+1):(High+Low),i])==0){
    immune[(High+1),i]=0.001
  }
  lowData=immune[1:High,i]
  highData=immune[(High+1):(High+Low),i]
  vioplot(lowData,at=3*(i-1),lty=1,add = T,col = '#E64B35')
  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col = '#00A087')
  wilcoxTest=wilcox.test(lowData,highData)
  p=wilcoxTest$p.value
  if(p<pFilter){
    cellPvalue=cbind(Cell=colnames(immune)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(lowData,highData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("High", "Low"),
       lwd=5,bty="n",cex=1.5,
       col=c("#E64B35","#00A087"))
text(seq(1,70,3),-0.1,xpd = NA,labels=colnames(immune),cex = 1,srt = 45,pos=2)
dev.off()

