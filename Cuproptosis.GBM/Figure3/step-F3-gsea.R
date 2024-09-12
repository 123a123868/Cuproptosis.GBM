# 画gsea的图：
# gase

# 
rm(list = ls())
options(stringsAsFactors = F)
library(clusterProfiler)
library(msigdbr)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(plyr)
library(dplyr)


## 数据处理函数
merge_result2 <- function(enrichResultList, output = "compareClusterResult") {
  if ( !is(enrichResultList, "list")) {
    stop("input should be a name list...")
  }
  if ( is.null(names(enrichResultList))) {
    stop("input should be a name list...")
  }
  x <- lapply(enrichResultList, as.data.frame)
  names(x) <- names(enrichResultList)
  y <- ldply(x, "rbind")   
  
  
  if (output == "compareClusterResult") {
    y <- plyr::rename(y, c(.id="Cluster"))
    y$Cluster = factor(y$Cluster, levels=names(enrichResultList))
    return(new("compareClusterResult",
               compareClusterResult = y))
  } 
  
  y <- plyr::rename(y, c(.id="Category"))
  if (output == "enrichResult") {
    return(new("enrichResult",
               result = y))        
  }
  
  if (output == "gseaResult") {
    return(new("gseaResult",
               result = y))        
  }   
}

keep_category <- function(em_ORA, n) {
  table_em <- as.numeric(table(em_ORA$Category))
  start <- rep(0, length(table_em) - 1)
  for(i in seq_len(length(table_em) - 1)) {
    start[i] <- sum(table_em[seq_len(i)])   
  }
  showCategorys <- sapply(table_em, function(x) min(n, x))
  start <- c(0, start) + 1
  end <- start + showCategorys - 1
  keep <- NULL
  for(i in seq_len(length(start))) {
    keep <- c(keep, c(start[i] : end[i]))
  } 
  return(keep)
}


enrich_filter <- function(em_result, showCategory) {
  keep <- keep_category(em_result, showCategory)
  em_result <- em_result[keep, ]
  if ("NES" %in% colnames(em_result))
    em_result$Count <- em_result$core_enrichment %>% 
    strsplit(split = "/")  %>%  
    vapply(length, FUN.VALUE = 1)
  return(em_result)
}

## 作图函数
em_plot <- function(em_1 = NULL, em_2 = NULL, showCategory = 2, fill = "p.adjust", hjust = 1) {
  
  fill <- match.arg(fill, c("Category", "p.adjust", "log10_p.adjust"))
  result1 <- enrich_filter(em_1, showCategory)     
  if (is.null(em_2)) { 
    result <- result1 
  } else {
    result2 <- enrich_filter(em_2, showCategory) 
    result2$Count <- -result2$Count          
    result <- rbind(result1, result2)
  }
  result$Category <- gsub("\n.*", "", result$Category)     
  result$log10_p.adjust <- log10(result$p.adjust)
  
  data_plot <- result[, c("ID", "Category", "p.adjust", "log10_p.adjust", "Count")]
  data_plot2 <- data_plot
  data_plot2$ID <- factor(data_plot2$ID, levels = unique(data_plot2$ID))
  data_plot2 <- plyr::rename(data_plot2, c("Count" = "gene_number"))
  h_just <- ifelse(data_plot2$gene_number < 0, -hjust, hjust)
  ggplot(data_plot2, aes_string(x = "gene_number", y = "ID", fill = fill)) + 
    geom_col() +       
    geom_text(aes_(x =~ gene_number + h_just, label =~ abs(gene_number)), 
              color="black") + 
    scale_x_continuous(label = abs,
                       expand = expansion(mult = c(.01, .01))) + #两侧留空
    theme_classic() + 
    ylab("") + 
    theme(axis.title.x = element_text(size = 15)) +     
    facet_grid(Category ~ ., scales="free", space="free") 
}

# 使用R包msigdbr里的注释
msigdbr_species()
gmt <- msigdbr(species = "Homo sapiens") 
gmt2 <- gmt%>%
  dplyr::select(gs_name, entrez_gene)
gmts <- split(gmt2, gmt$gs_cat)

# 加载差异表达分析结果
load("Figure3/out_data/deg.CC.Rdata")
# 区分上调和下调基因
## 设置阈值
logFC_t=0.58
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)
head(deg)
gsym.fc <- deg
gsym.fc$SYMBOL = rownames(gsym.fc)
gsym.fc[1:3,]
# 筛选差异基因
geneup = gsym.fc[gsym.fc$g == 'UP','SYMBOL'] 
genedown = gsym.fc[gsym.fc$g == 'DOWN','SYMBOL'] 
# 数量
length(geneup)
length(genedown)
# ORA富集分析
# 自定义函数
enrich_func <- function(x, gene, readable = FALSE) {
  en_result <- enricher(gene, TERM2GENE = x)
  if (readable & nrow(en_result) > 0) 
    en_result <- setReadable(en_result, 'org.Hs.eg.db', #物种
                             'ENTREZID')
  
  return(en_result)
}
# 把gene symbol转换为ENTREZ ID
# 此处物种是人，其他物种的ID转换方法，请参考FigureYa52GOplot
geneup.id <- bitr(geneup, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
# 富集分析
em_ORAup <- setNames(lapply(gmts, enrich_func, geneup.id$ENTREZID, readable = TRUE), names(gmts)) %>%
  merge_result2(output = "enrichResult")
genedown.id <- bitr(genedown, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
em_ORAdown <- setNames(lapply(gmts, enrich_func, genedown.id$ENTREZID, readable = TRUE), names(gmts)) %>%
  merge_result2(output = "enrichResult")
genesupdown <- union(geneup.id$ENTREZID, genedown.id$ENTREZID)
em_ORA <- setNames(lapply(gmts, enrich_func, genesupdown, readable = TRUE), names(gmts)) %>%
  merge_result2(output = "enrichResult") 

# 绘图
# 上调表达基因的富集结果
em_plot(em_ORAup, showCategory = 5, fill = "Category", hjust = 10) + scale_fill_brewer(palette="Set1")
# 下调表达基因的富集结果
em_plot(em_ORAdown, showCategory = 5, fill = "Category", hjust = 10) + scale_fill_brewer(palette="Set1")

# 上调的画在右侧，下调的画在左侧
em_plot(em_ORAup, em_ORAdown, showCategory = 3, fill = "Category", hjust = 12) + scale_fill_brewer(palette="Paired")
ggsave("Figure3/out_polt//top3_batchGSEA.pdf", width = 12, height = 15)

# 上调的画在右侧，下调的画在左侧
em_plot(em_ORAup, em_ORAdown, showCategory = 5, fill = "Category", hjust = 10) + scale_fill_brewer(palette="Paired")
ggsave("Figure3/out_polt//top5_batchGSEA.pdf", width = 12, height = 18)
