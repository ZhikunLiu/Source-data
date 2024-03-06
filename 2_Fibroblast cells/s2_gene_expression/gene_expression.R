library(Seurat)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)

# 读取数据
scRNA = read_rds("final_fibroblast cells.rds")

#增殖,凋亡,侵袭(迁移),胶原合成,TGF-β/Smad 信号通路

# 增殖：epithelial cell proliferation
# 凋亡：Bcl-2, Bax, Caspase 3
# 迁移：epithelial cell migration
# 胶原：FGF2, COL-I, COL-III
# TGF-β/Smad 信号通路：TGF-β,p-Smad 2,Smad 2,p-Smad 3,Smad3

# 增殖
scRNA@meta.data$new_group = paste0(scRNA@meta.data$cell_type,scRNA@meta.data$group)

bp = read.csv("bp.csv", header = T,, row.names = 1)
head(bp)

proliferation = unlist(strsplit(subset(bp, Description == "epithelial cell proliferation")$geneID, "/"))
pro_gene = bitr(proliferation, fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
head(pro_gene)


colors = c("#4ed5c7","#f7a67b","#d44865","#b8da8d")
p1 = DoHeatmap(scRNA, features = pro_gene$SYMBOL, group.by = "new_group",group.colors=colors,label=FALSE) +
  labs(title = "Epithelial cell proliferation") +
  scale_fill_gradient2(low = "#ffffe5",mid = "#d9f0a3",high = "#78c679")
p1


# 凋亡
p2 = DoHeatmap(scRNA, features = c('BCL2','BAX'), group.by = "new_group",group.colors=colors,label=FALSE) +
  labs(title = "Epithelial cell proliferation") +
  scale_fill_gradient2(low = "#ffffe5",mid = "#d9f0a3",high = "#78c679")

# 迁移
migration = unlist(strsplit(subset(bp, Description == "epithelial cell migration")$geneID, "/"))
pro_gene = bitr(migration, fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
head(pro_gene)
p3 = DoHeatmap(scRNA, features = pro_gene$SYMBOL, group.by = "new_group",group.colors=colors,label=FALSE) +
  labs(title = "Epithelial cell migration") +
  scale_fill_gradient2(low = "#ffffe5",mid = "#d9f0a3",high = "#78c679")
p3


# TGFβ信号通路:response to transforming growth factor beta
tgf = unlist(strsplit(subset(bp, Description == "response to transforming growth factor beta")$geneID, "/"))
pro_gene = bitr(tgf, fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
head(pro_gene)
p4 = DoHeatmap(scRNA, features = pro_gene$SYMBOL, group.by = "new_group",group.colors=colors,label=FALSE) +
  labs(title = "Response to transforming growth factor beta") +
  scale_fill_gradient2(low = "#ffffe5",mid = "#d9f0a3",high = "#78c679")
p4

plot_grid(p1,p2,p3,p4, ncol = 2)
