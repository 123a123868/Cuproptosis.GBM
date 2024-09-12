# 8-黑色素瘤  PRJEB23709
# 75名患者  9例PD1_EDT   41例PD1_PRE  9例ipiPD1_EDT  32例ipiPD1_PRE

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

# 读入数据
load("Figure7/8_weixin_PRJEB23709/PRJEB23709_Melanoma.Rdata")

# 检查数据：
dim(count.data) # 91个样本，符合
count.data[1:4,1:4]  #基因id已经转换，还没log
dim(phe)
head(phe)
# 用的是 Overall.Survival.Days 的生存时间

# 处理数据
# 标准化矩阵
# expr=log2(edgeR::cpm(count.data)+1)#标准化矩阵
# expr[1:4,1:4]


expr = log2(count.data+1)
expr[1:4,1:4]

# library(limma)
# expr = normalizeBetweenArrays(expr)

expr = as.data.frame(expr)
boxplot(expr$PD1_8_EDT)

# 临床信息
colnames(expr)
colnames(phe)
table(phe$sample_ID %in% colnames(expr))
rownames(phe) = phe$sample_ID
phe = phe[colnames(expr),]
identical(rownames(phe),colnames(expr))
write.csv(phe,file = "F7_drug/other_cohort/8_weixin_PRJEB23709/weixin_PRJEB23709.csv")

# 取32例联合治疗的PRE组
table(phe$group)
phe = phe[grep("ipiPD1_PRE",phe$group),]
dim(phe)
colnames(phe)
phe = phe[,c("Best.RECIST.response",
                 "Progression.Free.Survival.Days",
                 "Overall.Survival.Days",
                 "Last.Followup.Status")]
colnames(phe) = c("response","PFS","OS","event")
head(phe)
phe$response_2 = ifelse(phe$response =="PD","SD/PD","CR/PR")

phe$PFS = phe$PFS/30  # 换为月
phe$OS = phe$OS/30  # 换为月
boxplot(phe$OS)

table(phe$event)
phe$event_2 = ifelse( phe$event  %in% "Alive","0","1")
class(phe$event_2)
phe$event_2 = as.numeric(phe$event_2)


# 筛选矩阵
expr_2 = expr[,rownames(phe)]
dim(expr_2)
expr_2[1:4,1:4]

# 保存数据
save(phe,expr_2,file = "Figure7//8_weixin_PRJEB23709/ipiPD1_PRE_count.Rdata")


#######################################################################################################

#######################################################################################################
# 画图
# 设置
rm(list = ls())
load("Figure7//8_weixin_PRJEB23709/ipiPD1_PRE_count.Rdata")
expr_2[1:4,1:4]
head(phe)


#铁死亡基因
#cu死亡基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
table(cu_gene %in% rownames(expr_2))
expr_cu = expr_2[cu_gene,]

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
identical(rownames(data.score),rownames(phe))
data_all = cbind(data.score,phe)
res.cut <- surv_cutpoint(data_all, #数据集
                         time = "OS", #生存状态
                         event = "event_2", #生存时间
                         variables = "cu_score" #需要计算的数据列名
)
summary(res.cut) #查看数据最佳截断点及统计量
group <- as.factor(ifelse(data_all$cu_score > 0.7149917,'High','Low'))
data_all$group <- group 
table(group)

####  生存分析
sfit <- survfit(Surv(PFS, event_2)~group, data=data_all)
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

pdf('Figure7/8_weixin_PRJEB23709//8_a_survival_PFS.pdf' , onefile = F)
print(survp)
dev.off() 
        # 
        # # 多个变量画图
        # sfit1 <- survfit(Surv(PFS, event_2)~group, data=data_all)
        # sfit2 <- survfit(Surv(OS, event_2)~group, data=data_all)
        # fitlist <- list(sfit1, sfit2)
        # 
        # # 用ggsurvplot_combine画到一张图中
        # library(RColorBrewer)
        # color = brewer.pal(4,'Set1')  # 多种颜色
        # survp=ggsurvplot_combine(
        #   fitlist,    
        #   data = data_all,  # survfit object with calculated statistics.
        #   legend.title = 'group', 
        #   # legend = "top",#图例位置
        #   # legend.labs = c('High', 'Low'),
        #   pval = T, #在图上添加log rank检验的p值
        #   risk.table = TRUE, 
        #   risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
        #   xlab = "Time in Month", #x轴标题
        #   # xlim = c(0, 10), #展示x轴的范围
        #   break.time.by = 5, #x轴间隔
        #   size = 1.5, #线条大小
        #   ggtheme = theme_ggstatsplot(),
        #   palette= color, #配色
        # )
        # print(survp)
        # ggsurvplot_combine(fitlist, data = lung)
        # 


# 
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
table(data_all$response_2)

a <- data.frame(table(data_all$group,data_all$response_2))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
a$Freq

pvalue <- chisq.test(c(5, 2,21,4,ncol=2))$p.value #卡方检验
round(pvalue,digits = 2)
library(plyr)
ggplot(a,aes(Var1,percent,fill=Var2))+
  geom_bar(stat="identity",position = position_stack())+
  scale_fill_manual(values = c(brewer.pal(3,'Set1')),label=c("CR/PR","SD/PD"))+
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
ggsave(filename = "Figure7/8_weixin_PRJEB23709//8_CR_PR_Percent.pdf")

### 用小提琴图画
# 绘图数据
data = data_all
table(data$response_2)
data$cu_score
# 设置颜色
jco <- c("#008ECA","#DB423E") 
ggplot(data = data,aes(x = response_2, y = cu_score, fill = response_2))+
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
                     comparisons = list(c('CR/PR',"SD/PD")
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
ggsave("Figure7/8_weixin_PRJEB23709//8_c.pdf", width = 3, height = 4)

### 画ROC曲线
library(pROC)
library(ROCR)

# 绘图数据
data_all$cu_score
data_all$response_2

pred <- prediction(data_all$cu_score, data_all$response_2) 
#ROCR.simple$predictions为预测标签，ROCR.simple$labels为真实标签
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("fe_score (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

pdf("Figure7/8_weixin_PRJEB23709//8_ROC.pdf")
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("fe_score (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()
