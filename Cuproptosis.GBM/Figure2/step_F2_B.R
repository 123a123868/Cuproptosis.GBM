# maftool 突变分析


#设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)
library(maftools)       

# 读入10铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene

#
pdf(file="Figure2/out_polt/F2_B.pdf", width=8, height=6)
maf=read.maf(maf="Input_data/TCGA-GBM_MuTect2/input.maf")
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
print(vc_cols)
oncoplot(maf=maf, color=vc_cols, genes=cu_gene, draw_titv=T)
dev.off()

