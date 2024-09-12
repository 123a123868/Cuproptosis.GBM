# I．	切割小提琴图展示high和low两组间MHC分子表达差异: HLA-A, HLA-B, HLA-C, HLA-DMA,
# HLA-DMB, HLA-DOA, HLA-DOB, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DQB2, 
# HLA-DRA, HLA-DRB6, HLA-E, HLA-F, HLA-F-AS1, HLA-G, HLA-J.

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(ggplot2)

# 读入铁死亡评分high和low两组
load("Figure6/out_data/cu_score_group.Rdata")
head(new_dat)



# 读入表达矩阵
load("Figure3/out_data/expr.gbm.Rdata")
dim(expr.gbm)


# 匹配合并一下：
HLA_gene =  c("HLA-A","HLA-AS1","HLA-B","HLA-C","HLA-DMA","HLA-DMB" ,
              "HLA-DOA","HLA-DOB","HLA-DPA1",
              "HLA-DPB1","HLA-DPB2","HLA-DQA1","HLA-DQA2","HLA-DQB1",
              "HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-E",    
              "HLA-F","HLA-F-AS1","HLA-H","HLA-J","HLA-K","HLA-L","HLA-V",
              "HLA-DQB","HLA-DRB6","HLA-G") 

# 基因太多，展示不好看，删除部分没有差异的基因
HLA_gene =  c("HLA-A","HLA-C",
              "HLA-DOB","HLA-DQA2",
              "HLA-E",    
              "HLA-F","HLA-F-AS1",
              "HLA-DRB6","HLA-G") 

HLA_gene[!HLA_gene %in% rownames(expr.gbm)]


# [1] "HLA-AS1" "HLA-H"   "HLA-K"   "HLA-DQB" 没找到
data_HLA = expr.gbm[HLA_gene,]
data_HLA = na.omit(data_HLA)
dim(data_HLA)
data_HLA = as.data.frame(t(data_HLA))
new_dat = new_dat[rownames(data_HLA),]
identical(rownames(data_HLA),rownames(new_dat))
data_HLA$cu_group = new_dat$cu_group

#转换长数据
data_new = melt(data_HLA)
colnames(data_new) = c("group","gene","expression")


# 开始画小提琴图
#5.加载绘图函数
source("Figure6//Function_for_violin_plot.R")

#6.绘制小提琴图
## 6. 绘图
# 6.1 这里注意到原图用的是误差线，这里用步骤三加载的函数，计算一下误差信息
library(Rmisc)
Data_summary <- summarySE(
  data_new, measurevar="expression", groupvars=c("group","gene"))
head(Data_summary)


# 4.2. 出图
# 这个是我自己写的一个ggplot2的主题，可以自定义修改其中的参数
# 4.2. 出图
# 这个是我自己写的一个ggplot2的主题，可以自定义修改其中的参数
if(T){
  mytheme <- theme(plot.title = element_text(size = 8,color="black",hjust = 0.5),
                   axis.title = element_text(size = 8,color ="black"), 
                   axis.text = element_text(size= 8,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 8),
                   legend.title= element_text(size= 8)
  ) 
}

# 自行调整下面的参数
gene_split_violin <- ggplot(data_new,aes(x= gene,y= expression,fill= group))+
  geom_split_violin(trim= F,color="white",scale = "area") + #绘制分半的小提琴图
  geom_point(data = Data_summary,aes(x= gene, y= expression),pch= 80,
             position=position_dodge(0.5),size= 0.9)+ #绘制均值为点图
  geom_errorbar(data = Data_summary,aes(ymin = expression-ci, ymax= expression+ci), 
                width= 0.05, 
                position= position_dodge(0.5), 
                color="black",
                alpha = 0.8,
                size= 0.1) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ 
  labs(y=("Log2 expression"),x=NULL,title = "Split violin") + 
  theme_bw()+ mytheme +
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "anova",
                     label.y = max(data_new$expression),
                     hide.ns = T)
gene_split_violin;
ggsave(gene_split_violin,
       filename = "01-CC/violin_HLA-2.pdf",height = 6,width = 8)



# 4.2. 出图
# 这个是我自己写的一个ggplot2的主题，可以自定义修改其中的参数
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) 
}

# 自行调整下面的参数
gene_split_violin <- ggplot(data_new,aes(x= gene,y= expression,fill= group))+
  geom_split_violin(trim= F,color="white",scale = "area") + #绘制分半的小提琴图
  geom_point(data = Data_summary,aes(x= gene, y= expression),pch=19,
             position=position_dodge(0.5),size= 1)+ #绘制均值为点图
  geom_errorbar(data = Data_summary,aes(ymin = expression-ci, ymax= expression+ci), 
                width= 0.05, 
                position= position_dodge(0.5), 
                color="black",
                alpha = 0.8,
                size= 0.5) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ 
  labs(y=("Log2 expression"),x=NULL,title = "Split violin") + 
  theme_bw()+ mytheme +
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "anova",
                     label.y = max(data_new$expression),
                     hide.ns = T)
gene_split_violin;
ggsave(gene_split_violin,
       filename = "Figure6/out_polt//violin_HLA.pdf",height = 6,width = 8)
