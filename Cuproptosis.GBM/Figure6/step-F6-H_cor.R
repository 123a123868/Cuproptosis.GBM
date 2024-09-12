# 然后H图我想放的是这种，就是你不是算了
#Immune score、stroma score、ESTIMATE score、purity、TIDE、Dysfunction、Exclusion和TIS这8个score吗？想这个图这样，中间放的是FerrScore，然后两边各放4个score，计算FerrScore与另外8个评分之间的相关性，
#而不是直接去比较high和low两组之间各种评分的差异。

# 设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
library(limma)
library(survival)
library(ggplot2)


# 读入数据

# 读入铁死亡评分high和low两组
load("Figure6/out_data/cu_score_group.Rdata")
head(new_dat)

# 倒入免疫评分数据
load("Figure4/out_data/estimate_scores.Rdata")

# 匹配一下：
new_dat = new_dat[rownames(data_all),]
identical(rownames(data_all),rownames(new_dat))
data_all$cu_score = new_dat$cu_score
colnames(data_all)
data_all = data_all[,c("StromalScore","ImmuneScore","ESTIMATEScore","purity","cu_score",
                       "TIDE","Dysfunction","Exclusion","TIS_mean")]

# 求相关性，画热图
library(corrplot)
cor<- cor(data_all)
class(cor)
cor[1:5,1:5]
corrplot(cor, method = "square")
col2 <- colorRampPalette(c("#FFFF00","white" ,"#FF69B4"),alpha = TRUE)
corrplot(cor, addgrid.col = "grey70",type = "upper",
         outline = "orange",col = col2(100),method = "circle",diag = F,
         tl.col = "black")

#添加显著性标记：
#使用cor.mtest做显著性检验；
res1 <- cor.mtest(data_all, conf.level = .95)
res2 <- cor.mtest(data_all, conf.level = .99)

#提取p值矩阵；
p.mat = res1$p
p.mat[1:5,1:5]
corrplot(cor,col = col2(100),method = "color",
)

corrplot(cor, addgrid.col = "grey70",type = "upper",
         col = col2(100),method = "circle",diag = F,
         tl.col="black",tl.cex = 0.8,cl.pos = "r",cl.ratio = 0.2,
         p.mat = res1$p, sig.level = c(.001, .01, .05),outline="orange",
         insig = "label_sig",pch.cex = 1.2, pch.col = "black")


# 保存图片
pdf("Figure6/out_polt/H_cor.pdf",width = 10,height = 10)
values = c("#DB423E","#008ECA")
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                           "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                           "#D6604D", "#B2182B", "#67001F"))
corrplot(cor, addgrid.col = "grey70",type = "upper",
         col = col3(200),method = "circle",diag = F,
         tl.col="black",tl.cex = 0.8,cl.pos = "r",cl.ratio = 0.2,
         p.mat = res1$p, sig.level = c(.001, .01, .05),outline="orange",
         insig = "label_sig",pch.cex = 1.2, pch.col = "black")

dev.off()



# 保存数据
write.csv(cor,file = "F6_lasso/polt_lasso/H_cor.csv")
options(digits = 2)
options(scipen = 1000)
#使用cor.mtest做显著性检验；
res1 <- cor.mtest(data_all, conf.level = .95)
res2 <- cor.mtest(data_all, conf.level = .99)

#提取p值矩阵；
p.mat = res1$p
p.mat[1:5,1:5]

write.csv(p.mat,file = "F6_lasso/polt_lasso/H_p_value.csv")
