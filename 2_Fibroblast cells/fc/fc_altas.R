# 加载包
library(tidyverse)
library(Seurat)
library(patchwork)
library(harmony)
library(cols4all)
library(glmGamPoi) # SCT标准化加速包
library(cowplot)
library(clustree)



# 读取数据
scRNA = read_rds("final_scRNA.rds")
table(scRNA@meta.data$cell_type)

# 提取成纤维细胞
fc = subset(scRNA, cell_type == "Fibroblast cells")
fc@meta.data$cell_type = droplevels(fc@meta.data$cell_type)
# saveRDS(fc, "fibroblast cells.rds")

scRNA = readRDS("fibroblast cells.rds")


# SCT标准化可以直接跳过三步法标准化:NormalizeData, FidnVariableFeatures, ScaleData 
scRNA = SCTransform(scRNA, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = F)
# saveRDS(scRNA, file = "./sct_scRNA.rds")


# 读取SCT标准化后的数据
scRNA = readRDS("sct_fibroblast cells.rds")


# PCA降维
scRNA <- RunPCA(scRNA)
colnames(scRNA@meta.data)[1] = c("patient_ID", "nCount_RNA", "nFeature_RNA", "percent.mt") #重命名
# 去除批次效应
scRNA <- RunHarmony(scRNA, group.by.vars="patient_ID", max.iter.harmony=20, lambda=0.5, assay.use = "SCT")
# 查看主成分拐点
ElbowPlot(scRNA,ndims = 50)


# 使用前15个主成分进行聚类分析
scRNA <- FindNeighbors(scRNA, dims=1:15, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:15, reduction="harmony")

# 画图展示批次效应去除的结果
DimPlot(scRNA, reduction="umap", group.by="patient_ID", pt.size=1)+
  theme(legend.position="right", plot.title=element_blank())

# 使用聚类树选择合适的分辨率
obj = FindClusters(scRNA, resolution = seq(0.1, 1,by=0.1))
p4 = clustree(obj)
p4

# 根据聚类结果，分辨率大小设置为0.1。resolution = 0.1
scRNA <- FindClusters(scRNA, resolution=0.1)
p5 = UMAPPlot(scRNA, pt.size=1, label=T)+NoLegend()
p5

# 计算差异基因，并保存。
# scRNA.markers <- FindAllMarkers(scRNA, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# write.table(scRNA.markers, "markers.txt", sep = "\t")



# 被分为4个cluster
# 注释
cell_label = c("Fibroblast cells 1","Fibroblast cells 2", "Fibroblast cells 1", "Fibroblast cells 1")
# 给细胞的标签命名 
names(cell_label) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, cell_label)
scRNA[["cell_type"]] = Idents(scRNA)
colors = c("#4ed5c7","#f7a67b","#d44865","#b8da8d")
p8 = UMAPPlot(scRNA, pt.size=1, label=T, label.size=5, cols = colors)+
  NoLegend()
p8


# 计算每个样本中的各组细胞占比
cell_count = scRNA@meta.data %>%
  group_by(group, cell_type ) %>%
  count() %>% 
  group_by(group) %>% 
  mutate(Percent=n/sum(n))
head(cell_count)

# write.csv(cell_count, "cellcount_sample.csv", row.names = F)

cell_count = read.csv("cellcount_sample.csv", header = T)
p11 = ggplot(cell_count, aes(reorder(cell_type, -Percent*100), Percent*100, fill = group)) +
  geom_bar(stat = "identity",position = position_dodge()) +
  theme_bw()+
  geom_text(aes(label = round(Percent*100,2)), vjust = -0.5, hjust = 0.5,position = position_dodge(width = 0.8))+
  scale_fill_manual(values = c("#ff7b00","#006fbf")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top") +
  labs(title = "CellRatio (%)", x = "",y = "CellRatio", fill = "Sample")
p11

# saveRDS(scRNA, "final_fibroblast cells.rds")

head(scRNA)












# 差异基因富集分析
library(clusterProfiler)
library(org.Hs.eg.db)

scRNA = read_rds("final_fibroblast cells.rds")
table(scRNA@meta.data$cell_type)
# 使用findmarkers寻找差异基因
deg = FindMarkers(scRNA, ident.1 = "Fibroblast cells 2", ident.2 = "Fibroblast cells 1")
deg = arrange(deg, -avg_log2FC)
head(deg)
tail(deg)
dim(deg)

deg_fib1 = filter(deg,avg_log2FC<0)
deg_fib2 = filter(deg,avg_log2FC>0)

# 转换基因id
gene_fib1 = bitr(rownames(deg_fib1), fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
gene_fib2 = bitr(rownames(deg_fib2), fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

bp_fib1 = enrichGO(gene_fib1$ENTREZID,
         OrgDb = "org.Hs.eg.db",
         keyType = "ENTREZID",
         ont = "BP",
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.2,
         pAdjustMethod = "BH")

bp_fib2 = enrichGO(gene_fib2$ENTREZID,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   pAdjustMethod = "BH")

# write.csv(bp_fib1, "bp_fib1.csv")
# write.csv(bp_fib2, "bp_fib2.csv")

bp_fib1 = read.csv("bp_fib1_select.csv", header = T)
bp_fib2 = read.csv("bp_fib2_select.csv", header = T)

p12 = ggplot(bp_fib1, aes(reorder(Description, Count), Count)) +
  geom_col(aes(fill = pvalue)) + 
  coord_flip() +
  scale_color_gradient(low = "#4ed5c7",high = "#d2d8ce")+
  scale_fill_gradient(low = "#4ed5c7",high = "#d2d8ce")+
  theme_bw() +
  geom_text(aes(y = 0, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = -0.03) + # 左对齐
  theme(axis.text.y = element_blank()) + #主题中去掉y轴通路标签
  labs(title = "Fibroblast cells 1", x = "Description")
p12


p13 = ggplot(bp_fib2, aes(reorder(Description, Count), Count)) +
  geom_col(aes(fill = pvalue)) + 
  coord_flip() +
  scale_color_gradient(low = "#f7a67b",high = "#d2d8ce")+
  scale_fill_gradient(low = "#f7a67b",high = "#d2d8ce")+
  theme_bw() +
  geom_text(aes(y = 0, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = -0.03) + # 左对齐
  theme(axis.text.y = element_blank()) + #主题中去掉y轴通路标签
  labs(title = "Fibroblast cells 2", x = "Description")
p13



library(cowplot)
plot_grid(p12,p13,ncol = 1)

tail(deg_fib1)
head(deg_fib2)
rownames(deg_fib1)
rownames(deg_fib2)

# fib1 genes: SFRP2, FBLN1, COL1A1, COL1A2, COL3A1, LUM
# fib2 genes: RGS5, ID4, MCAM, TINAGL1, MT1A, PDGFA

genes = c('SFRP2', 'FBLN1', 'COL1A1', 'COL1A2', 'COL3A1', 'LUM',
          'RGS5', 'ID4', 'MCAM', 'TINAGL1', 'MT1A', 'PDGFA')
DotPlot(scRNA, features = genes, group.by = "cell_type", cols = c("#006fbf","#ff7b00")) +
  theme_bw() +
  coord_flip()


