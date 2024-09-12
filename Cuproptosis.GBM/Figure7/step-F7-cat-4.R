
#  anti-CTLA4

# 设置
rm(list = ls()) 
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)
library(ggplot2)
library(GEOquery)
library(AnnoProbe)
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
library(clusterProfiler)


# 读入数据
expr_CAT4 = read.csv("Figure7/tide_cat4//Nathanson2017_CTLA4_Melanoma_RNASeq_Post/ICB.Nathanson2017_Ipilimumab_Melanoma_Post.self_subtract",
                     header = T,fill = T,sep = '\t')
phe = read.csv("Figure7/tide_cat4//Nathanson2017_CTLA4_Melanoma_RNASeq_Post/ICB.Nathanson2017_Ipilimumab_Melanoma_Post.clinical",
               header = T,fill = T,sep = '\t')

# 检查数据
boxplot(expr_CAT4$P4631)  # 经过sacle的数据
# 转换id
gene.id <- bitr(expr_CAT4$Entrez, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
row.names(expr_CAT4) = expr_CAT4$Entrez
expr_CAT4 = expr_CAT4[gene.id$ENTREZID,]
identical(rownames(expr_CAT4),gene.id$ENTREZID)
row.names(expr_CAT4) = gene.id$SYMBOL
expr_CAT4 = expr_CAT4[,-1]
expr_CAT4[1:4,1:4]


#铜死亡基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
table(cu_gene %in% rownames(expr_CAT4))

# 临床信息：
phe$Response = ifelse(phe$Response =="1","Response","no_Response")
table(phe$Response)
boxplot(phe$OS) #年的时间
table(phe$OS.Event)


# 保存数据
save(phe,expr_CAT4,file = "Figure7/tide_cat4//count.Rdata")


# 画图
# 设置
rm(list = ls())
load("Figure7/tide_cat4//count.Rdata")

#铜死亡基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
expr_cu = expr_CAT4[cu_gene,]

data.score = as.data.frame(t(expr_cu))
data.score$cu_score = exp(-0.262490103*data.score$CDKN2A
                          -0.266355127*data.score$DLAT
                          +0.184258117*data.score$DLD
                          -0.374608653*data.score$FDX1
                          -0.086511329*data.score$GLS
                          -0.242718907*data.score$LIAS
                          +0.508961778*data.score$LIPT1
                          +0.220024283*data.score$MTF1
                          +0.001332578*data.score$PDHA1
                          +0.149181298*data.score$PDHB)
##    surv_cutpoint函数，花费高低分组的最佳切点
library(survminer)
rownames(phe) = phe$Patient
identical(rownames(data.score),rownames(phe))
data_all = cbind(data.score,phe)
colnames(data_all)
res.cut <- surv_cutpoint(data_all, #数据集
                         time = "OS", #生存状态
                         event = "OS.Event", #生存时间
                         variables = "cu_score" #需要计算的数据列名
)
summary(res.cut) #查看数据最佳截断点及统计量
group <- as.factor(ifelse(data_all$cu_score > 0.9608801,'High','Low'))
data_all$group <- group 
table(group)

####  生存分析
sfit <- survfit(Surv(OS, OS.Event)~group, data=data_all)
sfit
summary(sfit)
## more complicate figures.
survp=ggsurvplot(
  sfit,                     # survfit object with calculated statistics.
  legend.title = 'group', 
  # legend = "top",#图例位置
  # legend.labs = c('High', 'Low'),
  pval = T, #在图上添加log rank检验的p值
  risk.table = TRUE, 
  risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
  xlab = "Time in Month", #x轴标题
  # xlim = c(0, 10), #展示x轴的范围
  break.time.by = 5, #x轴间隔
  size = 1.5, #线条大小
  ggtheme = theme_ggstatsplot(),
  palette= c("#DB423E","#008ECA"), #配色
)
print(survp)

pdf('Figure7/tide_cat4//a_survival_os.pdf' , onefile = F)
print(survp)
dev.off() 

##################################
# 高低分组的，反应。
#######
#比例图：
library(dplyr)
library(plyr)
library(ggplot2)
library(RColorBrewer)

# 绘图数据
table(data_all$group)
table(data_all$Response)

a <- data.frame(table(data_all$group,data_all$Response))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
a$Freq

pvalue <- chisq.test(c(5, 3,6,1,ncol=2))$p.value #卡方检验
round(pvalue,digits = 2)
library(plyr)
ggplot(a,aes(Var1,percent,fill=Var2))+
  geom_bar(stat="identity",position = position_stack())+
  scale_fill_manual(values = c(brewer.pal(3,'Set1')),label=c("no_Response","Response"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="Cu-Score",y="Percent Weidght",
       fill="")+
  geom_text(aes(label=label),vjust=3,size=6,color="black")+
  annotate(geom = "text",
           cex=6,
           x=1.5, y=105, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(pvalue<0.001, "<0.001 ", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
ggsave(filename = "Figure7/tide_cat4//CR_PR_Percent.pdf")

### 用小提琴图画
# 绘图数据
data = data_all
table(data$Response)
data$cu_score
# 设置颜色
jco <- c("#008ECA","#DB423E") 
ggplot(data = data,aes(x = Response, y = cu_score, fill = Response))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),
              size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = FALSE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21, size=1, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("fe-score")) +
  xlab("group")  +
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = list(c('no_Response',"Response")
                     ),
                     method="wilcox.test"
  )+#添加检验
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))
# 保存图像
ggsave("Figure7/tide_cat4//8_c.pdf", width = 3, height = 4)

### 画ROC曲线
library(pROC)
library(ROCR)

# 绘图数据
data_all$cu_score
data_all$Response

pred <- prediction(data_all$cu_score, data_all$Response) 
#ROCR.simple$predictions为预测标签，ROCR.simple$labels为真实标签
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("cu_score (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

pdf("Figure7/tide_cat4//8_ROC.pdf")
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("cu_score (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()












