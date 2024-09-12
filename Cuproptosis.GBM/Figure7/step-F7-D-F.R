##  可视化IMvigor结果

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)
library(data.table)
library(dplyr)
library(futile.logger) 
library(glmSparseNet)
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)

# 读入数据
load('Figure7/out_data/IMvigor_count.Rdata')
load("Figure7/out_data/cu_scroe_group.Rdata")
names(pdata)[names(pdata) == 'time'] <- 'time.m'
names(pdata)[names(pdata) == 'event'] <- 'status'
identical(row.names(dat_cox),rownames(pdata))
dat_cox = cbind(dat_cox,pdata)

#  绘图
#  Immune_phenotype
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#作图
table(dat_cox$Immune_phenotype)
dat_cox2 = dplyr::filter(dat_cox, !is.na(Immune_phenotype))
ggplot(dat_cox2,aes(x=Immune_phenotype,y=cu_score,fill=Immune_phenotype))+
  geom_boxplot()+
  theme_bw()+#改变绘图主题+
  theme(
    panel.grid = element_blank(),#去掉背景网格
    legend.position = c('none')#去掉图例
  )+scale_fill_manual(values=c(brewer.pal(3,'Set1')))+
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = list(c('desert',"excluded"),
                                        c('desert',"inflamed"),
                                        c('excluded',"inflamed")),
                     method="wilcox.test"
  )+#添加检验
  xlab("Immune_phenotype")#修改横坐标
ggsave(filename = "Figure7/out_polt/F7-D.pdf")

###  Enrollment_IC
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#作图
table(dat_cox$Enrollment_IC)
ggplot(dat_cox,aes(x=Enrollment_IC,y=cu_score,fill=Enrollment_IC))+
  geom_boxplot()+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = c('none')
  )+scale_fill_manual(values=c(brewer.pal(3,'Set1')))+
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = list(c('IC0',"IC2"),
                                        c('IC0',"IC1"),
                                        c('IC1',"IC2")),
                     method="t.test"
  )+
  xlab("IC_Level")
ggsave(filename = "Figure7/out_polt/F7-E.pdf")

##  TC_Level 
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#作图
table(dat_cox$TC_Level)
dat_cox = dplyr::filter(dat_cox, !is.na(TC_Level))
ggplot(dat_cox,aes(x=TC_Level,y=cu_score,fill=TC_Level))+
  geom_boxplot()+
  theme_bw()+#改变绘图主题+
  theme(
    panel.grid = element_blank(),#去掉背景网格
    legend.position = c('none')#去掉图例
  )+scale_fill_manual(values=c(brewer.pal(3,'Set1')))+
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = list(c('TC0',"TC1"),
                                        c('TC0',"TC2+"),
                                        c('TC1',"TC2+")),
                     method="t.test"
  )+#添加检验
  xlab("TC_Level")#修改横坐标
ggsave(filename = "Figure7/out_polt/F7-F.pdf")

