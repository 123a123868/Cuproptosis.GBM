# 热图可视化
# 读入数据，合并
# 然后可视化

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
library(glmSparseNet)
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)

# 读入数据
load("Figure4/免疫浸润多算法/out_data/ssgsea_results.Rdata")
ssgseaOut = ssgseaOut[-1,]
rownames(ssgseaOut) 
rownames(ssgseaOut)  = paste0(rownames(ssgseaOut),"_ssGSEA")
rownames(ssgseaOut)


load("Figure4/免疫浸润多算法/out_data/ciber.results.Rdata")
colnames(ciber.results)
ciber.results = as.data.frame(ciber.results)
ciber.results = ciber.results[,!colnames(ciber.results) %in% c("P-value","Correlation","RMSE")]
ciber.results = as.data.frame(t(ciber.results))
rownames(ciber.results)
rownames(ciber.results)  = paste0(rownames(ciber.results),"_CIBERSORT")

# 匹配
ciber.results = ciber.results[,colnames(ssgseaOut)]
identical(colnames(ciber.results),colnames(ssgseaOut))



load("Figure4/免疫浸润多算法/out_data/mcp.results.Rdata")
rownames(mcp.results)
rownames(mcp.results)  = paste0(rownames(mcp.results),"_MCPCOUNTER")
identical(colnames(ciber.results),colnames(mcp.results))


### 
load("Figure4/免疫浸润多算法/out_data/xcell_results.Rdata")
rownames(xcell_results)
rownames(xcell_results)  = paste0(rownames(xcell_results),"_XCELL")
rownames(xcell_results)
identical(colnames(ciber.results),colnames(xcell_results))
xcell_results = xcell_results[,colnames(ciber.results)]
identical(colnames(ciber.results),colnames(xcell_results))

##
timer_results = read.csv("Figure4/免疫浸润多算法/out_data/timer_ruslt.csv")
rownames(timer_results) = timer_results$sampleID
timer_results = timer_results[,-1]
timer_results = t(timer_results)
rownames(timer_results)
rownames(timer_results) <- gsub('_',' ',rownames(timer_results))
rownames(timer_results)
rownames(timer_results)  = paste0(rownames(timer_results),"_TIMER")
rownames(timer_results)
identical(colnames(ciber.results),colnames(timer_results))
timer_results = timer_results[,colnames(ciber.results)]
identical(colnames(ciber.results),colnames(timer_results))

## 全部合并在一起
data_all = rbind(timer_results,
                 ciber.results,
                 ssgseaOut,
                 xcell_results,
                 mcp.results)


#######
# 开始画热图

# 读入高低风险的分组
# 读入铁死亡评分high和low两组
load("Figure6/out_data/cu_score_group.Rdata")
head(new_dat)

# 设计不同算法的注释
annotation_row = data.frame(
  cell = rownames(data_all),
  group = "na"
)
annotation_row$group = str_split(annotation_row$cell,'[_]',simplify = T)[,2]
rownames(annotation_row) = annotation_row$cell
annotation_row = data.frame(
  group2 = annotation_row$group,
  row.names = annotation_row$cell
)



##### 画热图
# 画图
write.csv(data_all,file = "Figure4/免疫浸润多算法/out_data/data_all.csv")
data_all = read.csv("Figure4/免疫浸润多算法/out_data/data_all.csv")
row.names(data_all) = data_all$X
data_all = data_all[,-1]


pheatmap(data_all,show_colnames =F,show_rownames = F) 
n=t(scale(t(data_all))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)

#排序一下
new_dat = new_dat[colnames(n),]
new_dat = new_dat[order(new_dat$cu_score),]
annotation_col = data.frame(group=new_dat$cu_group)
rownames(annotation_col) = rownames(new_dat)
annotation_col$group
n = n[,rownames(annotation_col)]
pheatmap(n,show_colnames =F,show_rownames = T,
         cluster_row = F,cluster_col = F,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         cellheight = 8.7,
         filename = 'Figure4/免疫浸润多算法/out_polt/heatmap_all.pdf')
dev.off()


cellheight = 8.7
pheatmap(n,show_colnames =F,show_rownames = T,
         annotation_col= F,
         annotation_row = annotation_row,
         cluster_row = F,cluster_col = F)




pheatmap(n,show_colnames =F,show_rownames = T,
         annotation_col=ac,
         annotation_row = annotation_row,
         cluster_row = F,cluster_col = F,
         cellheight = 8.7,
         filename = 'Figure4/out_polt/F4_A.pdf')
dev.off()





n=t(scale(t(data_all))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
#排序一下
pd.CC = pd.CC[order(pd.CC$cluster),]
ac=data.frame(group=pd.CC$cluster)
rownames(ac)=pd.CC$sample
ac$group
n = n[,rownames(ac)]
pheatmap(n,show_colnames =F,show_rownames = T,
         annotation_col=ac,
         annotation_row = annotation_row,
         cluster_row = F,cluster_col = F,
         cellheight = 8.7,
         filename = 'Figure4/out_polt/F4_A.pdf')
dev.off()
































