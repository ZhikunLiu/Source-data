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
scRNA = read_rds("scRNA.rds")
#向scRNA新增一列percent.mt数据
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

VlnPlot(
  scRNA, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
  # cols = colors,
  pt.size = 0.1, 
  ncol = 3
)


# 线粒体基因比例<10%, 细胞内检测基因个数200-6000
scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<6000&percent.mt<20)

VlnPlot(
  scRNA, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  # cols = colors,
  pt.size = 0, 
  ncol = 3
)


# SCT标准化可以直接跳过三步法标准化:NormalizeData, FidnVariableFeatures, ScaleData 
# scRNA = SCTransform(scRNA, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = F)
# saveRDS(scRNA, file = "./sct_scRNA.rds")

# 读取SCT标准化后的数据
scRNA = readRDS("sct_scRNA.rds")


# PCA降维
scRNA <- RunPCA(scRNA)
colnames(scRNA@meta.data)[1] = c("patient_ID", "nCount_RNA", "nFeature_RNA", "percent.mt") #重命名
# 去除批次效应
scRNA <- RunHarmony(scRNA, group.by.vars="patient_ID", max.iter.harmony=20, lambda=0.5, assay.use = "SCT")
# 查看主成分拐点
ElbowPlot(scRNA,ndims = 50)

# 使用前20个主成分进行聚类分析
scRNA <- FindNeighbors(scRNA, dims=1:20, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:20, reduction="harmony")

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
head(scRNA)
scRNA = subset(scRNA, seurat_clusters %in% c(0,1,2,3,4,5,6))
scRNA@meta.data$seurat_clusters = droplevels(scRNA@meta.data$seurat_clusters)
UMAPPlot(scRNA, pt.size=1, label=T)+NoLegend()
head(scRNA)

# 计算差异基因，并保存。
# scRNA.markers <- FindAllMarkers(scRNA, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# write.csv(scRNA.markers, "markers.csv")

# 注释
VlnPlot(scRNA, features=c("CPA3"), pt.size=0,group.by = "seurat_clusters")+
  NoLegend()+
  theme(axis.title.x=element_blank())

# cluster 0: Fibraoblast cells:THY1,TAGLN,MYH11,FGF7,ASPN,GPC3,SFRP1,WISP2,MFAP5
# cluster 1: Endothelial cells: "FLT1", "KDR", "VWF","PECAM1", "CDH5"
# cluster 2: Epidermal stem cells: KRT5, KRT14, TP63,KRT1,SBSN,KRTDAP
# cluster 3: fibroblast cells:
# cluster 4: Macrophages: CD163, CSF1R, CYBB, FPR3
# cluster 5: Endothelial cells:
# cluster 6: Mast cells: "SLC18A2","CPA3","HPGDS","TPSB2"


fc = c('THY1','TAGLN','MYH11','FGF7','ASPN','GPC3','SFRP1','WISP2','MFAP5')
ec = c("FLT1", "KDR","PECAM1", "CDH5","VWF")
esc = c('KRT5', 'KRT14', 'TP63','KRT1','SBSN','KRTDAP')
mac = c('CD163', 'CSF1R', 'CYBB', 'FPR3')
mc = c("SLC18A2","CPA3","HPGDS","TPSB2")
anno_markers = c(fc,ec,esc,mac,mc)

# cluster 注释
cell_label = c("Fibroblast cells","Endothelial cells","Epidermal stem cells",
               "Fibroblast cells","Macrophages","Endothelial cells","Mast cells")
# 给细胞的标签命名 
names(cell_label) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, cell_label)
scRNA[["cell_type"]] = Idents(scRNA)

# UMAP注释绘图
colors = c("#5fc9f8","#53d769","#fd9426","#fc3158","#147efb")
p8 = UMAPPlot(scRNA, pt.size=1, label=T, label.size=5, cols = colors)+
  NoLegend()
p8


# 绘制点图
p9 = DotPlot(scRNA, features=anno_markers, cols=c("#006fbf", "#ff7b00"))+coord_flip()+
  theme(
    axis.text.x=element_text(angle=30, hjust=1, size=10, face="bold"), 
    axis.text.y=element_text(face="bold", size=12), 
    axis.title.x=element_blank(), 
    axis.title.y=element_blank()
  )
p9


genes = c("MFAP5","WISP2","PECAM1","KRT5","CSF1R","CPA3")
# 绘制小提琴图
p10 = VlnPlot(scRNA, features=genes, pt.size=0, ncol = 3, cols = colors)+
  NoLegend()+
  theme(axis.title.x=element_blank())
p10


# 添加细胞分组
head(scRNA@meta.data)
scRNA@meta.data$group = ifelse(scRNA@meta.data$patient_ID %in% c("GSM4994379","GSM4994380","GSM4994381"),"keloid","normal scar")


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

# saveRDS(scRNA, "final_scRNA.rds")

