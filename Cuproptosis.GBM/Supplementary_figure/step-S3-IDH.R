# 8个基因在不同比较IDH_wildtype和IDH_mutant 10个铜基因的差异

# 设置
rm(list = ls())  
options(stringsAsFactors = F) 
library(stringr)
library(stringi)
library(data.table)

# 读入10个铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene

# 导入表达矩阵数据
#读入数据
load('Figure3/out_data/CGGA_count.Rdata')
dim(expr.gbm)
data_gene = expr.gbm[cu_gene,]
data_gene = as.data.frame(t(data_gene))

# 临床信息
phe = fread("Input_data/CGGA_GBM/CGGA.mRNAseq_693_clinical.20200506.txt/CGGA.mRNAseq_693_clinical.20200506.txt",
              header = T, sep = '\t',data.table = F)
phe[1:4,1:4]
rownames(phe) = phe$CGGA_ID
phe = phe[colnames(expr.gbm),]
identical(rownames(phe),rownames(data_gene))
colnames(phe)
data_gene$PRS_type = phe$PRS_type


#加载包
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#作图
colnames(data_gene)

#Grade
data_gene =  dplyr::filter(data_gene, !is.na(PRS_type))
table(data_gene$PRS_type)
data_gene = na.omit(data_gene)
table(data_gene$PRS_type)

for (i in 1:10) {
  # i = 4
  gene = cu_gene[i]
  data = data_gene[c(gene,"PRS_type")]
  colnames(data) = c("gene","Grade")
  ggplot(data,aes(x=Grade,y=gene,color=Grade))+
    geom_boxplot()+
    theme_bw()+#改变绘图主题+
    theme(
      panel.grid = element_blank(),#去掉背景网格
      legend.position = c('none')#去掉图例
    )+
    scale_color_manual(values=brewer.pal(3,'Set1'))+#修改颜色
    stat_compare_means(aes(label = ..p.signif..),comparisons = list(c('Primary','Recurrent')),method="wilcox.test"
    )+#添加检验
    ylab(gene)#修改横坐标
  ggsave(filename = paste0("Supplementary_figure_3/out_polt/PRS_type//Grade_",gene,".pdf"))
}

#######################################
##################
# T期的
colnames(data_all)
data_all_T =  dplyr::filter(data_all, !is.na(T))
table(data_all_T$T)

for (i in 1:8) {
  # i = 4
  gene = gene_fe[i]
  data = data_all_T[c(2,4+i)]
  colnames(data) = c("Tgrate","gene")
  data$Tgrate = paste0('stage',data$Tgrate) 
  ggplot(data,aes(x=Tgrate,y=gene,color=Tgrate))+
    geom_boxplot()+
    theme_bw()+#改变绘图主题+
    theme(
      panel.grid = element_blank(),#去掉背景网格
      legend.position = c('none')#去掉图例
    )+
    scale_color_manual(values=brewer.pal(3,'Set1'))+#修改颜色
    stat_compare_means(aes(label = ..p.signif..),comparisons = list(c('stage2','stage1'),
                                                                    c('stage3','stage2'),
                                                                    c('stage3','stage1')),method="t.test"
    )+#添加检验
    ylab(gene)#修改横坐标
  ggsave(filename = paste0("F_S3/T/T_",gene,".pdf"))
}


######################
# N期
colnames(data_all)
data_all_N =  dplyr::filter(data_all, !is.na(N))
table(data_all_N$N)
data_all_N = data_all_N[!data_all_N$N == "1mi",]
table(data_all_N$N)

for (i in 1:8) {
  # i = 4
  gene = gene_fe[i]
  data = data_all_N[c(3,4+i)]
  colnames(data) = c("Ngrate","gene")
  data$Ngrate = paste0('stage',data$Ngrate) 
  table(data$Ngrate)
  ggplot(data,aes(x=Ngrate,y=gene,color=Ngrate))+
    geom_boxplot()+
    theme_bw()+#改变绘图主题+
    theme(
      panel.grid = element_blank(),#去掉背景网格
      legend.position = c('none')#去掉图例
    )+
    scale_color_manual(values=brewer.pal(4,'Set1'))+#修改颜色
    stat_compare_means(aes(label = ..p.signif..),comparisons = list(c('stage0','stage1'),
                                                                    c('stage0','stage2'),
                                                                    c('stage0','stage3'),
                                                                    c('stage1','stage2'),
                                                                    c('stage1','stage3'),
                                                                    c('stage2','stage3')
    ),method="t.test"
    )+#添加检验
    ylab(gene)#修改横坐标
  ggsave(filename = paste0("F_S3/N/N_",gene,".pdf"))
}



