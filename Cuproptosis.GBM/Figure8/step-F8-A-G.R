# F8图片-A

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

# Responder
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#作图
colour = c("#DB423E","#7CFC00")
table(dat_cox$pd1)
ggplot(dat_cox,aes(x=pd1,y=cu_score,fill=pd1))+
  geom_boxplot()+
  theme_bw()+#改变绘图主题+
  theme(
    panel.grid = element_blank(),#去掉背景网格
    legend.position = c('none')#去掉图例
  )+scale_fill_manual(values= colour)+
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = list(c('NonResponder',"Responder")),
                     method="wilcox.test"
  )+#添加检验
  xlab("Anti-PDL1-Responder")#修改横坐标
ggsave(filename = "Figure8/out_polt/F8-A.pdf")

################################################################
# CR_PR_Percent
library(dplyr)
library(plyr)
table(dat_cox$Best_Confirmed_Overall_Response)
a <- data.frame(table(dat_cox$Best_Confirmed_Overall_Response,dat_cox$cu_group))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
a = a[!a$Freq == 0,] 

pvalue <- chisq.test(c(25,43,63,167,ncol=2))$p.value 
library(plyr)
ggplot(a,aes(Var1,percent,fill=Var2))+
  geom_bar(stat="identity",position = position_stack())+
  scale_fill_manual(values = c(brewer.pal(4,'Set1')),label=c("High","Low"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+ 
  labs(x="Type",y="Percent Weidght",
       fill="")+
  geom_text(aes(label=label),vjust=3,size=6,color="black")+
  annotate(geom = "text",
           cex=6,
           x=1.5, y=105, 
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), 
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
ggsave(filename = "Figure8/out_polt/F8-B.pdf")

#############################################################
# CR_PR_Percent
library(dplyr)
library(plyr)
a <- data.frame(table(dat_cox$cu_group,dat_cox$BOR_binary))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
a$Freq
pvalue <- chisq.test(c(13,91,55,139,ncol=2))$p.value 
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
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
ggsave(filename = "Figure8/out_polt/F8-C.pdf")





#############################################################
# TMB
# 用小提琴图画
data = dat_cox
table(data$cu_group)
colnames(data)
data$FMOne_mutation_burden_per_MB
data = dplyr::filter(data, !is.na(FMOne_mutation_burden_per_MB))
# 设置颜色
jco <- c("#008ECA","#DB423E") 
ggplot(data = data,aes(x = cu_group, y = FMOne_mutation_burden_per_MB, fill = cu_group))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),
              size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = FALSE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21, size=1, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("FMOne_mutation_burden_per_MB")) +
  xlab("group")  +
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = list(c('High',"Low")
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
ggsave("Figure8/out_polt//8_D.pdf", width = 3, height = 4)




#############################################################
# TNB
# 用小提琴图画
data = dat_cox
table(data$cu_group)
colnames(data)
data$Neoantigen_burden_per_MB
data = dplyr::filter(data, !is.na(Neoantigen_burden_per_MB))
# 设置颜色
jco <- c("#008ECA","#DB423E") 
ggplot(data = data,aes(x = cu_group, y = Neoantigen_burden_per_MB, fill = cu_group))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),
              size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = FALSE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21, size=1, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("Neoantigen_burden_per_MB")) +
  xlab("group")  +
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = list(c('High',"Low")
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
ggsave("Figure8/out_polt//8_E.pdf", width = 3, height = 4)



#########################
# km
sfit <- survfit(Surv(time, event)~cu_group, data=dat_cox)
sfit
summary(sfit)
## more complicate figures.
survp=ggsurvplot(
  sfit,                     # survfit object with calculated statistics.
  legend.title = 'Risk level', 
  # legend = "top",#图例位置
  # legend.labs = c('High', 'Low'),
  pval = T, #在图上添加log rank检验的p值
  risk.table = TRUE, 
  risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
  xlab = "Time in years", #x轴标题
  # xlim = c(0, 10), #展示x轴的范围
  break.time.by = 1, #x轴间隔
  size = 1.5, #线条大小
  ggtheme = theme_ggstatsplot(),
  palette= c("#DB423E","#008ECA"), #配色
)
print(survp)
pdf('Figure8/out_polt/F8-G.pdf' , onefile = F)
print(survp)
dev.off() 


### 画ROC曲线
library(pROC)
library(ROCR)

##########################
#  绘图数据
dat_cox$cu_score
dat_cox$BOR_binary

pred <- prediction(dat_cox$cu_score, dat_cox$BOR_binary) 
#ROCR.simple$predictions为预测标签，ROCR.simple$labels为真实标签
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("cu_score (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

pdf("Figure8/out_polt//8_ROC.pdf")
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("cu_score (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()



