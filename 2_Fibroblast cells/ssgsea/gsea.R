# devtools::install_github("junjunlab/GseaVis")
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
library(tidyverse)
library(dplyr)
library(GseaVis)


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



# 读取单细胞数据
scRNA = readRDS("final_fibroblast cells.rds")

# 使用findmarkers寻找差异基因
deg = FindMarkers(scRNA, ident.1 = "Fibroblast cells 1", ident.2 = "Fibroblast cells 2")
deg = arrange(deg, -avg_log2FC)
# write.table(deg,"deg.txt", sep = "\t")
deg = read.table("deg.txt", header = T, row.names = 1, sep = "\t")
head(deg)

#symbol转entrez ID：
symbol <- rownames(deg)
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)

#准备genelist文件(entrez+log2FC)：
genelist <- deg$avg_log2FC
names(genelist) <- rownames(deg)
#过滤掉ID转换中缺失部分基因：
genelist <- genelist[names(genelist) %in% entrez[,1]]
head(genelist)
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
head(genelist)

#按照log2FC从高到低排序：
genelist <- sort(genelist,decreasing = T)
head(genelist)


#GSEA_KEGG富集分析：
R.utils::setOption( "clusterProfiler.download.method","auto") ##如果富集时报错就加上这句代码
KEGG_ges <- gseGO(
  geneList = genelist,
  ont = "BP",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
#将entrez重新转换为symbol：
KEGG_ges2 <- setReadable(KEGG_ges,
                         OrgDb = org.Hs.eg.db,
                         keyType="ENTREZID")
#转换前：
KEGG_ges@result$core_enrichment[1]
#转换后：
KEGG_ges2@result$core_enrichment[1]


save(KEGG_ges,KEGG_ges2,file = c('GSEA_KEGG.Rdata'))
rm(list=ls())
########
#2.使用GseaVis包完成GSEA结果可视化：
load('GSEA_KEGG.Rdata')
result <- KEGG_ges2@result
View(result) #查看GSEA结果表，选择想展示的基因集/通路绘图

# 保存文件
# write.csv(result,"result.csv")

colnames(result)

dim(result)
#2.1——标准GSEA富集分析结果图：
#标注Top rank基因：
#2.1——标准GSEA富集分析结果图：
gseaNb(
  object = KEGG_ges2,
  geneSetID = KEGG_ges2@result$ID[88], #绘制结果表中第38个pathway(自行选择即可)
  subPlot = 3, #常规为3图组合，如果不需要条码图或rank图可以设置为1 or 2
  addPval = T, #图中是否显示P值和NES标签
  pvalX = 0.95,
  pvalY = 0.8 #调整标签X/Y坐标控制位置
)


#更改配色：
gseaNb(
  object = KEGG_ges2,
  geneSetID = KEGG_ges2@result$ID[88],
  subPlot = 3,
  addPval = T,
  pvalX = 0.95,
  pvalY = 0.8,
  curveCol = c('#7582c1', '#dd568d'), #ES折线图颜色更改
  htCol = c("#7582c1", "#dd568d"), #热图条颜色更改
  rankCol = c("#7582c1", "white", "#dd568d") #rank分布图颜色更改
)


#去掉散点和热图条：
gseaNb(
  object = KEGG_ges2,
  geneSetID = KEGG_ges2@result$ID[88],
  addPval = T,
  pvalX = 0.95,
  pvalY = 0.8,
  newGsea = T,
  addPoint = F, #是否添加散点
  rmHt = T #是否移除底部热图条
)


#更改配色：
gseaNb(
  object = KEGG_ges2,
  geneSetID = KEGG_ges2@result$ID[88],
  addPval = T,
  pvalX = 0.95,
  pvalY = 0.8,
  newGsea = T,
  addPoint = F,
  rmHt = T, #是否移除底部热图条
  newCurveCol = c("#4ed5c7","#b9b4ad", "#f7a67b"), #ES点阵图颜色更改
  newHtCol = c("#4ed5c7", "white", "#f7a67b"), #热图条颜色更改
  addGene = T, #是否添加基因
  markTopgene = T, #是否标注Top基因
  topGeneN = 20, #标注前多少个gene
  kegg = T, #是否将entrez转symbol
  geneCol = '#4d4d4d', #基因名标签颜色更改
  rmSegment = T #是否移除红线
)



