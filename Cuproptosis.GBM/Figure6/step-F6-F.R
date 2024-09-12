# 桑葚图

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)

# 读入数据
load("Figure3/out_data/expr.gbm.Rdata")
load("Figure3/out_data/pd.CC.Rdata")
dim(expr.gbm)
identical(colnames(expr.gbm),pd.CC$sample)
load("Figure6/out_data/cu_score_group.Rdata")
head(new_dat)

# 合并
new_dat$sample = rownames(new_dat)
new_dat = new_dat[pd.CC$sample,]
all.data = cbind(new_dat,pd.CC)
all.data$event = ifelse(all.data$event == 1,"Death","Alive")
colnames(all.data)
all.data = all.data[,c("event","cu_group","cluster")]

# 画图
library(ggalluvial)
library(RColorBrewer)
library(ggplot2)
library(ggalluvial)
gg <- ggplot(all.data,
             aes(axis1 = cluster, axis2 = cu_group,axis3 = event ))
gg
gg <- gg +
  geom_alluvium(aes(fill = as.factor(cluster)),
                width = 2/5, discern = FALSE);gg

gg <- gg +geom_stratum(width = 2/5, discern = FALSE);gg
#添加文本
gg <- gg +geom_text(stat = "stratum", discern = FALSE,
                    aes(label = after_stat(stratum))) ;gg
gg <- gg +
  theme(legend.position = "none",#去除刻度线和背景颜色
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size =15,colour = "black"),#坐标轴名
        axis.title = element_blank()) +
  scale_x_discrete(position = "top");gg #坐标轴位置
ggsave('Figure6/out_polt/F6_F.pdf')
