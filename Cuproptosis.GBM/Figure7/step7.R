# 7-黑色素瘤
# GSE78220
# GPL11154
# anti-PD-1  28个病人

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

# 下载数据：
# 读入数据
library(openxlsx)
expr<- read.xlsx("Figure7/7/GSE78220_PatientFPKM.xlsx",sheet=1)
expr[1:4,1:4]
rownames(expr) = expr$Gene
expr = expr[,-1]
expr[1:4,1:4]
boxplot(expr$Pt1.baseline)
# 标准化矩阵
expr = log2(expr+1)
boxplot(expr$Pt1.baseline)
# 修改列名
library(stringr)
colnames(expr) <- gsub('.baseline','',colnames(expr))
colnames(expr)
expr[1:4,1:4]

#铜死亡基因
# cu评分
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
table(cu_gene %in% rownames(expr))




# 临床信息
    # gset <- getGEO('GSE78220', destdir=".",
                   AnnotGPL = T,     ## 注释文件
                   getGPL = T)       ## 平台文件
    ggset <- gset[[1]]
    pdata <- pData(ggset) 
    write.csv(pdata,file = "F7_drug/other_cohort/7/7_pdata.csv")
pdata = read.csv("Figure7/7/7_pdata.csv")
colnames(pdata) 
pdata = pdata[,c("title","characteristics_ch1.6","vital.status.ch1","anti.pd.1.response.ch1")]
head(pdata)
rownames(pdata) = pdata$title
pdata$title[!pdata$title %in% colnames(expr)]
colnames(expr)[14] = "Pt16"
table(pdata$title %in% colnames(expr))
pdata$time = substr(pdata$characteristics_ch1.6,start=26,stop=30)
pdata = pdata[!pdata$time =="NA",]
pdata$time = as.numeric(pdata$time)
pdata$time = pdata$time/30  # 换为月
boxplot(pdata$time)
pdata$event = ifelse( pdata$vital.status.ch1 %in% "Dead","1","0")
class(pdata$event)
pdata$event = as.numeric(pdata$event)

pdata$anti.pd.1.response.ch1
a = substr(x= str_split_fixed(string = pdata$anti.pd.1.response.ch1,pattern = ' ',n=2)[,1],start=1,stop=1)
b = substr(x= str_split_fixed(string = pdata$anti.pd.1.response.ch1,pattern = ' ',n=2)[,2],start=1,stop=1)
paste0(a,b)
pdata$treament = paste0(a,b)

pdata = pdata[,c(5,6,7)]
head(pdata)

save(pdata,expr,file = "Figure7//7/count.Rdata")



#####################################################################################
# 设置
rm(list = ls())
load("Figure7/7/count.Rdata")
expr[1:4,1:4]
head(pdata)

# z-score表达谱
#cexpr <- t(scale(t(expr))) 
expr[1:4,1:4]
dim(expr)


#cu死亡基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
table(cu_gene %in% rownames(expr))
expr_cu = expr[cu_gene,]

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
identical(rownames(data.score),rownames(pdata))
data.score = data.score[rownames(pdata),]
identical(rownames(data.score),rownames(pdata))
data.score$time = pdata$time
data.score$event = pdata$event


class(pdata$event)
res.cut <- surv_cutpoint(data.score, #数据集
                         time = "time", #生存状态
                         event = "event", #生存时间
                         variables = "cu_score" #需要计算的数据列名
)
summary(res.cut) #查看数据最佳截断点及统计量
group <- as.factor(ifelse(data.score$fe_score > 1.325757,'High','Low'))
data.score$group <- group 
table(group)

####  生存分析
sfit <- survfit(Surv(time, event)~group, data=data.score)
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

pdf('Figure7/7//7_a_survival.pdf' , onefile = F)
print(survp)
dev.off() 



# 
##################################
# 高低分组的，反应。
#######
#比例图：
library(dplyr)
library(plyr)
library(ggplot2)
library(RColorBrewer)
identical(rownames(data.score),rownames(pdata))
data.score$treament = pdata$treament
data.score$pd1_treament = ifelse(data.score$treament =="PD","SD/PD","CR/PR")

a <- data.frame(table(data.score$group,data.score$pd1_treament))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
a$Freq

pvalue <- chisq.test(c(1, 3,13,10,ncol=2))$p.value #卡方检验
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
           label=paste0("P ", ifelse(pvalue<0.001, " < 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
ggsave(filename = "Figure7/7/7_CR_PR_Percent.pdf")


### 用小提琴图画
# 绘图数据
data = data.score
table(data$pd1_treament)
data$fe_score
# 设置颜色
jco <- c("#008ECA","#DB423E") 
ggplot(data = data,aes(x = pd1_treament, y = fe_score, fill = pd1_treament))+
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
ggsave("Figure7/7/7_c.pdf", width = 3, height = 4)



### 画ROC曲线
library(pROC)
library(ROCR)
data(ROCR.simple)
pred <- prediction(data.score$fe_score, data.score$pd1_treament) 
#ROCR.simple$predictions为预测标签，ROCR.simple$labels为真实标签
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("fe_score (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

pdf("Figure7/7/7_ROC.pdf")
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("fe_score (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()


