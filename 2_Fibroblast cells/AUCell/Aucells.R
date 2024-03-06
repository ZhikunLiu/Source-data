# 加载包
library(tidyverse)
library(Seurat)
library(cols4all)
library(glmGamPoi) # SCT标准化加速包
library(cowplot)
library(clustree)
library(AUCell)
library(clusterProfiler) # 读取gmt文件
library(ggsignif)
library(GSVA)
library(ggpubr)
library(org.Hs.eg.db)


# 读取文件
scRNA = readRDS("final_fibroblast cells.rds")
table(scRNA@meta.data$cell_type)


#Build gene expression rankings for each cell
CellRank<- AUCell_buildRankings(as.matrix(scRNA@assays$RNA@data))


# TGF-βsignaling pathway:hsa04350
# 下载信号通路基因集
library(KEGGREST)
# 获取基因列表
kegg_gene = keggGet("hsa04350")

gene = sapply(strsplit(kegg_gene[[1]]$GENE, " "), "[", 1)

data_gene = data.frame("Term" = rep(kegg_gene[[1]]$PATHWAY_MAP, length(gene)),
                       "Gene" = gene)
head(data_gene)
gene_id = bitr(gene, fromType = "ENTREZID",toType = "SYMBOL", OrgDb = "org.Hs.eg.db", drop = T)

gene_id$term = rep(kegg_gene[[1]]$PATHWAY_MAP, length(gene_id$SYMBOL))
head(gene_id)
# write.csv(gene_id, "tgfb signaling pathway.csv", row.names = F)


signaling = split(gene_id$SYMBOL,data_gene$Term)
signaling

#Calculates the 'AUC' for each gene-set in each cell.
cells_AUC <- AUCell_calcAUC(signaling, CellRank,nCores = 5, aucMaxRank=nrow(CellRank)*0.05)
aucs <- getAUC(cells_AUC)[names(signaling), ]

# aucs %>%
#   as.data.frame() %>%
#   write.csv("aucs.csv")

aucs = read.csv("aucs.csv", header = T)
head(aucs)
names(aucs) = c("samples", "Score")
head(aucs)

head(scRNA@meta.data)
group = dplyr::select(scRNA@meta.data, c("group","cell_type"))
head(group)
group$samples = rownames(group)

# 合并数据
data = merge(aucs, group, by = "samples")
head(data)


table(data$group)

library(gghalves)
mycol = c("#ff7b00","#006fbf")
p1 = ggplot(subset(data, cell_type == "Fibroblast cells 1"), aes(reorder(cell_type, -Score), Score, fill = group)) +
  geom_half_point(aes(color = group),alpha = 0.9, side = 'L', size = 0.2) +
  geom_half_violin(side = 'R', alpha = 0.6) +
  geom_half_boxplot(alpha = 0.3, outlier.size = 0.1) +
  theme_classic() +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = c("grey","grey")) +
  coord_flip() +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5),legend.position = "top") +
  labs(title = "Fibroblast cells 1",x = "", y = "AUCells Score") +
  stat_compare_means(aes(group = group),label = "p.format",size=3,method = "t.test")
p1


p2 = ggplot(subset(data, cell_type == "Fibroblast cells 2"), aes(reorder(cell_type, -Score), Score, fill = group)) +
  geom_half_point(aes(color = group),alpha = 0.9, side = 'L', size = 0.2) +
  geom_half_violin(side = 'R', alpha = 0.6) +
  geom_half_boxplot(alpha = 0.3, outlier.size = 0.1) +
  theme_classic() +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = c("grey","grey")) +
  coord_flip() +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5),legend.position = "top") +
  labs(title = "Fibroblast cells 2",x = "", y = "AUCells Score") +
  stat_compare_means(aes(group = group),label = "p.format",size=3,method = "t.test")
p2


