# LASSO的图

# 设置
rm(list = ls())  
options(stringsAsFactors = F)
options(digits = 10)
library(data.table)
library(dplyr)
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


# 0. load data----
load('Figure6/out_data//coefs.v_lasso_model.Rdata')
xdata[1:4, 1:3]
# xdata = xdata[,-3]   #保留两个基因 
head(ydata)
# 1. multivarible-cox----
dat <- ydata[, -3]
head(dat)  
cg <- sort(names(coefs.v))
cg <- names(coefs.v)
#cg = cg[-3]
n <- apply(xdata,2,scale)[,cg] %>%
  as.data.frame()
head(n)
rownames(n) <- rownames(xdata)
boxplot(n)
dat_cox <- cbind(dat,n)
head(dat_cox)
multivariate <- paste(sort(colnames(n)), collapse = '+') 
attach(dat_cox)
s <-  paste0(' Surv(time, event) ~  ', multivariate )
model <- coxph(as.formula(s), data = dat_cox )
summary(model, data = dat_cox)

# 计算cu评分。用e的评分公式
# coefficient
# CDKN2A -0.157630 
# LIAS   -0.009714
data.score = as.data.frame(xdata)
model
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
identical(rownames(dat_cox),rownames(data.score))
dat_cox$cu_score = data.score$cu_score


# 2.  visualization----
## 2.1 forest----
library(survminer)
options(scipen=1)
ggforest(model, data = dat_cox, 
         main = "Hazard ratio", 
         cpositions = c(0.06, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)
ggsave('Figure6/out_polt/F6_H.pdf', height = 5, width = 13)


## ROC
library(timeROC)
head(dat_cox)
new_dat <- dat_cox[, c('event', 'time','cu_score')]
head(new_dat)
## need 3 cols，time、event and risk scores
result <- with(new_dat, timeROC(T=time,
                                delta=event,
                                marker=cu_score,
                                cause=1,
                                times = c(1, 3, 5),
                                iid = TRUE))
#identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))

## 将geom_line()改为geom_smooth(method = "loess")
## 在数学中有种算法叫“样条插补法”，这种方法可以获得过点的平滑曲线
ggplot() + 
  geom_smooth(data = dat, 
              aes(x = fpr, y = tpr, color = time), 
              size = 1,
              method = "loess",
              se = FALSE) + 
  scale_color_manual(name = NULL,
                     values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,3,5),"-year survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x = c(0, 1), y = c(0,1)), 
            color = "grey",
            linetype = 'dotdash')+
  theme_ggstatsplot()+
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.line = element_line(linetype = 1),
        panel.grid = element_blank(),
        legend.background = element_rect(linetype = 2, size = 0.2, colour = "black"),
        legend.position = c(0.665,0.135))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
ggsave('Figure6/out_polt/F6_G.pdf')



##  km----
## surv_cutpoint函数，高低分组的最佳切点
res.cut <- surv_cutpoint(new_dat, 
                         time = "time", 
                         event = "event", 
                         variables = "cu_score" 
)
summary(res.cut) #查看数据最佳截断点及统计量
cu_group <- as.factor(ifelse(new_dat$cu_score > 0.7012854027,'High','Low'))
new_dat$cu_group <- cu_group 
sfit <- survfit(Surv(time, event)~cu_group, data=new_dat)
sfit
summary(sfit)
save(new_dat,file ='Figure6/out_data/cu_score_group.Rdata')

## KM
survp=ggsurvplot(
  sfit,                 
  legend.title = 'cu_score', 
  legend.labs = c('High', 'Low'),
  pval = T, 
  risk.table = TRUE, 
  risk.table.y.text = F,
  xlab = "Time in years", 
  break.time.by = 1, 
  size = 1.5, 
  ggtheme = theme_ggstatsplot(),
  palette="nejm", 
)
print(survp)
pdf('Figure6/out_polt/F6_E.pdf', onefile = F)
print(survp)
dev.off()

## risk_socre
library(ggrisk)
library(rms)
dc <- datadist(dat_cox)
options(datadist="dc")
pdf('Figure6/out_polt/F6_B_D.pdf', onefile = F)
ggrisk(model,
       color.A = c(low = "#0B5E9D", high = "#EA5205"),
       color.B = c(code.0 = "#0B5E9D", code.1 = "#EA5205"),
       color.C = c(low = "#0B5E9D", median = "white", high = "#EA5205"))
dev.off() 
save(model,dat_cox,new_dat, file = 'Figure6/out_data//multicox_model.Rdata')



