## Setting options
rm(list=ls())
options(stringsAsFactors = F)


library(tidyverse)
library(clusterProfiler) # Y叔的包有没有，这个其实只是为了ID转换
library(org.Hs.eg.db) 
library(STRINGdb)
library(igraph)
library(ggraph)



## Import 10 cuproptosis-regulated genes
cu_gene = read.table("../Input_data//cu_gene.txt")
cu_gene = cu_gene$V1
cu_gene










## Load R packages
library(stringr)
library(igraph)


## Load interaction data from string (v11.5)
links= read.delim("string_interactions.tsv")
network <- graph_from_data_frame(d=links[,c(1:2,10)], directed=F) 
deg <- degree(network, mode="all")

nodes <- data.frame(
  name=unique(links$X.node1),
  group=c( rep("up",6),rep("down",2)))
network <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
my_color = c("#66C2A5", "#FC8D62", "#8DA0CB")[as.numeric(as.factor(V(network)$group))]
par(bg="grey13", mar=c(0,0,0,0))

plot(network,vertex.size=deg,
     layout=layout.circle,
     vertex.color=my_color,
     vertex.label.cex=0.7,
     vertex.label.color="white",
     vertex.frame.color="transparent",
     edge.width=E(network)$combined_score*3,
     edge.curved=0.1)
pdf('A2.pdf',width = 4, height = 4)

plot(network, 
     vertex.size=deg,
     layout=layout.ci,
     vertex.color=my_color,
     vertex.label.cex=0.7,
     vertex.label.color="white",
     vertex.frame.color="transparent",
     edge.width=E(network)$combined_score*3,
     edge.curved=0.1)

pdf('A_2.pdf',width = 4, height = 4)
plot(network, 
     vertex.size=deg,
     layout=layout.circle,
     vertex.color=my_color,
     vertex.label.cex=0.7,
     vertex.label.color="white",
     vertex.frame.color="transparent",
     edge.width=E(network)$combined_score*3,
     edge.curved=0.1)
dev.off()


