# E．	桑葚图，内容包括cluster，铁死亡评分(high/low)，Intrinsic_Subtype, 
# mRNA_Subtype, 
# 生存结局(dead/ alive)；


#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(ggalluvial)
library(RColorBrewer)
library(ggplot2)
library(ggalluvial)


# 读入CC分型信息
pd.CC = fread( 'Figure3/CC-result/CC-result.k=2.consensusClass.csv',
               header = F, data.table = F)
colnames(pd.CC) = c('sample','cluster')
pd.CC$cluster = paste0('C',pd.CC$cluster) 

# cu死亡评分
load('Figure6/out_data/cu_score_group.Rdata')
data_fe_score = new_dat
data_fe_score$event = ifelse(data_fe_score$event == 1,"Death","Alive")


# 读入临床信息
library(data.table)
phe = fread("Input_data/CGGA_GBM/CGGA.mRNAseq_693_clinical.20200506.txt/CGGA.mRNAseq_693_clinical.20200506.txt",
              header = T, sep = '\t',data.table = F)
phe[1:4,1:4]
rownames(phe) = phe$CGGA_ID
colnames(phe)
phe = phe[,c("PRS_type","IDH_mutation_status")]


# 合并称为一个大数据集
all.data = data_fe_score
rownames(pd.CC) = pd.CC$sample
all.data = all.data[rownames(pd.CC),]
identical(rownames(pd.CC),rownames(all.data))
all.data$CC_cluster = pd.CC$cluster

phe = phe[rownames(all.data),]
identical(rownames(phe),rownames(all.data))
all.data = cbind(all.data,phe)
write.csv(all.data,file = "Figure6/out_data/cu_phe_cc.csv")


# 转换长数据
colnames(all.data)
all.data = all.data[,c("CC_cluster","cu_group","event","PRS_type","IDH_mutation_status")]
UCB_lodes <- to_lodes_form(all.data[,1:ncol(all.data)],
                           axes = 1:ncol(all.data),
                           id = "Cohort")
dim(UCB_lodes)
head(UCB_lodes)
tail(UCB_lodes)


# 开始画桑基图
#画图
library(ggalluvial)
library(RColorBrewer)
library(ggplot2)
library(ggalluvial)
colnames(all.data)

library(ggalluvial)
#定义足够多的颜色，后面从这里选颜色
mycol <- rep(c("#088247","#D20A13","#FFD121",
               "#223D6C","#11AA4D","#58CDD9","#7A142C",
               "#5D90BA","#029149","#431A3D","#91612D","#6E568C",
               "#E0367A","#D8D155","#64495D","#7CC767"),2)
ggplot(UCB_lodes,
       aes(x = x, stratum = stratum, alluvium = Cohort,
           fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) + 
  geom_flow(width = 1/8) + #线跟方块间空隙的宽窄
  geom_stratum(alpha = .9,width = 1/4) + #方块的透明度、宽度
  geom_text(stat = "stratum", size = 4,color="black") + #文字大小、颜色
  
  #不喜欢默认的配色方案，用前面自己写的配色方案
  scale_fill_manual(values = mycol) +
  
  xlab("") + ylab("") +
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + #去掉坐标轴
  ggtitle("")+
  guides(fill = FALSE) 
ggsave("Figure6/out_polt/G_sanky.pdf")

