# ssgase
# 细胞免疫浸润

# 设置
rm(list = ls())
options(stringsAsFactors = F)
library(GSVA)
library(limma)
library(GSEABase)

# 读入数据
load("Figure3/out_data/expr.gbm.Rdata")
load("Figure3/out_data/pd.CC.Rdata")
dim(expr.gbm)
expr.gbm[1:4,1:4]
identical(colnames(expr.gbm),pd.CC$sample)

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
ssgseaOut[1:4,1:4]

# save data----
save(ssgseaOut,file = "Figure4/免疫浸润多算法/out_data/ssgsea_results.Rdata")
write.table(ssgseaOut,file="Figure4/out_data/ssgseaOut.txt",
            sep="\t",quote=F,col.names=F)
