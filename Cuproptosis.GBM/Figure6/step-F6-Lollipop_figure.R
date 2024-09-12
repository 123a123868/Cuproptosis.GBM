# F．	棒棒糖图展示8个基因与铁死亡评分的相关性；

#设置
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
library(ggplot2)

# 读入fe-socre评分和8基因的表达矩阵

# 铜死亡评分
load('Figure6/out_data/cu_score_group.Rdata')
index <- "cu_score" #基因名

# 计算相关性
# 读入10铜基因
cu_gene = read.table("Input_data/cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene
load("Figure3/out_data/expr.gbm.Rdata")
expr_cu = expr.gbm[cu_gene,]
expr_cu = as.data.frame(t(expr_cu))

# 合并两个数据集
new_dat = new_dat[rownames(expr_cu),]
identical(rownames(expr_cu),rownames(new_dat))
data_cu_score = expr_cu
data_cu_score$cu_score = new_dat$cu_score

# 计算相关性
y <- as.numeric(data_cu_score$cu_score)
head(y)

colnames <- colnames(data_cu_score)[1:10]
colnames
data <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(data_cu_score[,i]),y,type="spearman")
  data[i,2] <- test$estimate                                            
  data[i,3] <- test$p.value
}
names(data) <- c("symbol","correlation","pvalue")
head(data)

data %>% 
  #filter(pvalue <0.05) %>% # 如果不想把p值大于0.05的放在图上，去掉最前面的#号
  ggplot(aes(correlation,forcats::fct_reorder(symbol,correlation))) +
  geom_segment(aes(xend=0,yend=symbol)) +
  geom_point(aes(col=pvalue,size=abs(correlation))) +
  scale_colour_gradientn(colours=c("#DB423E","#008ECA")) +
  #scale_color_viridis_c(begin = 0.5, end = 1) +
  scale_size_continuous(range =c(2,8))  +
  theme_minimal() +
  ylab(NULL)

ggsave('Figure6/out_polt/F_cor_8gene_fescore_2.pdf')



