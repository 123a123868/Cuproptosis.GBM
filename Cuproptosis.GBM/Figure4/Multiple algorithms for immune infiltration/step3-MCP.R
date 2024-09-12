# MCP算法
# 免疫浸润
# 设置
rm(list=ls())
library(MCPcounter)
library(stringr)
library(data.table)
library(pheatmap)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# 读入数据
# TPM数据 。TPM 是RNAseq测序结果里很好的归一化表达矩阵 有正数和负数
# 直接log的表达矩阵就是可以的了。
# 读入数据
load("Figure3/out_data/expr.gbm.Rdata")
expr.gbm[1:4,1:4]
targetCancerTPMestimates <- MCPcounter.estimate(expr.gbm, featuresType = "HUGO_symbols")
dim(targetCancerTPMestimates)

# 保存数据
mcp.results = targetCancerTPMestimates
save(mcp.results,file = "Figure4/免疫浸润多算法/out_data/mcp.results.Rdata")
