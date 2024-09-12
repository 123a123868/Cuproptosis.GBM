# go 和 kegg  分析

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
options(digits = 4)

# 读入数据
load("Figure3/out_data/expr.gbm.Rdata")
load("Figure3/out_data/pd.CC.Rdata")
dim(expr.gbm)
identical(colnames(expr.gbm),pd.CC$sample)

# 标准化
library(limma)
expr.gbm=normalizeBetweenArrays(expr.gbm)
boxplot(expr.gbm[,1:6],las=2)

# 差异分析
dat = expr.gbm 
group_list = pd.CC$cluster
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
deg = topTable(fit,coef=2,adjust='BH', n=Inf) 
save(deg,file = "Figure3/out_data/deg.CC.Rdata")

# 火山图
library(EnhancedVolcano)
EnhancedVolcano(deg,
                lab =  rownames(deg),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05, 
                FCcutoff = 0.58,
                pointSize = 1.0, 
                labSize = 3.0
                )
ggsave(filename = 'Figure3/out_polt/EnhancedVolcano_for_DEG_DEseq2.pdf')


# 画图
# kegg
rm(list=ls())
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

# 读入数据
load("Figure3/out_data/deg.CC.Rdata")
head(deg)
## 设置阈值
logFC_t=0.58
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)
head(deg)
deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
head(DEG)
gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )
data(geneList, package="DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)

# 画图
# go和kegg

# 上调基因
gene_up=unique(gene_up)
gene_up
# go分析
go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all") 
library(ggplot2)
library(stringr)
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave("Figure3/out_polt/F3_F.pdf") 

# kegg分析
kk.up <- enrichKEGG(gene         = gene_up,
                    organism     = 'hsa',
                    #universe     = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk.up)[,1:6]
kk=kk.up
g_kegg = dotplot(kk)
g_kegg
ggsave(g_kegg,filename = "Figure3/out_polt/F3_H.pdf")

#########################################################################
# 下调基因
gene_down=unique(gene_down)
gene_down
# go分析
go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all") 
library(ggplot2)
library(stringr)
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave("Figure3/out_polt/F3_G.pdf") 

# kegg分析
kk.up <- enrichKEGG(gene         = gene_down,
                    organism     = 'hsa',
                    #universe     = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk.up)[,1:6]
kk=kk.up
g_kegg = dotplot(kk)
g_kegg
ggsave(g_kegg,filename = "Figure3/out_polt/F3_I.pdf")
