library(Seurat)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)


scRNA = read_rds("final_fibroblast cells.rds")

matrix = scRNA@assays$SCT@data %>% 
  as.data.frame()

# 读取bp
bp = read.csv("bp.csv", header = T)
head(bp)

gene = c(bp$GOBP_FIBROBLAST_PROLIFERATION[bp$GOBP_FIBROBLAST_PROLIFERATION!=""],
         bp$GOBP_FIBROBLAST_MIGRATION[bp$GOBP_FIBROBLAST_PROLIFERATION!=""],
         bp$GOBP_EXECUTION_PHASE_OF_APOPTOSIS[bp$GOBP_EXECUTION_PHASE_OF_APOPTOSIS!=""],
         bp$GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA[bp$GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA!=""],"TXNDC5")
length(gene) # 588

data = subset(matrix, rownames(matrix) %in% gene)
dim(data) # 443 20866

# 保存数据
# data %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   write.csv("gene_expr.csv")


data = read.csv("gene_expr.csv", header = T, row.names = 1)
head(data)
data$samples = rownames(data)

group = dplyr::select(scRNA@meta.data, c(cell_type, group))
group$samples = rownames(group)

data1 = merge(data, group, by = "samples")

data_long = pivot_longer(data1, cols = colnames(data)[!colnames(data) %in% c("TXNDC5","samples")], names_to = "gene",values_to = "expression")
head(data_long)

# 过滤0值
data2 = filter(data_long,TXNDC5 != 0 & expression !=0)
head(data2)
dim(data2)
dim(data_long)


mycol = c("#4ed5c7","#f7a67b")
p1 = ggplot(data2, aes(cell_type, TXNDC5))+
  geom_point(position = 'jitter', color = 'grey',
             size = 2, alpha = 0.8) +
  geom_violin(aes(fill = cell_type),
              color = 'grey', alpha = 0.8,
              scale = 'width',
              linewidth = 0.6, #外轮廓粗细
              trim = TRUE) +
  geom_boxplot(color = 'white',
               outlier.color = 'black',
               width = 0.4, #箱子宽度
               size = 0.8, #外轮廓描边粗细
               fill = NA) +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol) +
  theme_classic() +
  theme(legend.position = "top")
p1

mycol = c("#ff7b00","#006fbf")
p2 = ggplot(data2, aes(group, TXNDC5))+
  geom_point(position = 'jitter', color = 'grey',
             size = 2, alpha = 0.8) +
  geom_violin(aes(fill = group),
              color = 'grey', alpha = 0.8,
              scale = 'width',
              linewidth = 0.6, #外轮廓粗细
              trim = TRUE) +
  geom_boxplot(color = 'white',
               outlier.color = 'black',
               width = 0.4, #箱子宽度
               size = 0.8, #外轮廓描边粗细
               fill = NA) +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol) +
  theme_classic()+
  theme(legend.position = "top")
p2

plot_grid(p1,p2, ncol = 2, labels = c("A","B"))


# 标准化矩阵
head(data1)
data1 = subset(data1, cell_type == "Fibroblast cells 1")
data3 = dplyr::select(data1,-c(samples,cell_type,group))
head(data3)
expr_mean = colMeans(data3)
head(expr_mean)
expr_var = apply(data3,2,sd)
head(expr_var)
expr = (data3-expr_mean)/expr_var
head(data3)
expr$samples = data1[,"samples"]
expr = dplyr::select(expr, samples, everything())
head(expr)
# write.csv(expr, "expr(Normalized).csv", row.names = F)












# 计算AUC得分
library(clustree)
library(AUCell)
library(clusterProfiler) # 读取gmt文件
library(ggsignif)
library(GSVA)

# 查看bp
head(bp)
head(scRNA)

fc1 = subset(scRNA, cell_type == "Fibroblast cells 1")
fc1@meta.data$cell_type =droplevels(fc1@meta.data$cell_type)
#Build gene expression rankings for each cell
CellRank<- AUCell_buildRankings(as.matrix(fc1@assays$RNA@data))
bp_list = as.list(bp)
head(bp_list)
# 使用lapply结合匿名函数移除空字符串
bp_list <- lapply(bp_list, function(x) x[x != ""])
head(bp_list)

#Calculates the 'AUC' for each gene-set in each cell.
cells_AUC <- AUCell_calcAUC(bp_list, CellRank,nCores = 5, aucMaxRank=nrow(CellRank)*0.05)

# 要测试的基因集
geneSet <- unique(names(bp_list))
geneSet
aucs <- getAUC(cells_AUC)[geneSet, ]
head(aucs)

# t(aucs) %>%
#     as.data.frame() %>%
#     write.csv("aucs.csv")

aucs = read.csv("aucs.csv", header = T, row.names = 1)
head(aucs)
dim(aucs) # 14714     4
summary(aucs)
# 去除含有0的行
aucs_filtered <- aucs[!apply(aucs == 0, 1, any), ]
dim(aucs_filtered) # 14482     4
summary(aucs_filtered)
# write.csv(aucs_filtered,"aucs_filtered.csv")




# 计算相关性
library(correlation)
head(data1)
df = data1 %>% 
  dplyr::select(c(samples, cell_type, group, TXNDC5)) %>% 
  subset(cell_type == "Fibroblast cells 1")
head(df)
dim(df)
df_mean = mean(df$TXNDC5)
df_var = var(df$TXNDC5)
df_txndc5 = (df$TXNDC5-df_mean)/df_var
summary(df_txndc5)


aucs = read.csv("aucs.csv", header = T, row.names = 1)
head(aucs)
aucs_select = subset(aucs, rownames(aucs) %in% df$samples)
dim(aucs_select)

cor_df = correlation(df_txndc5, aucs_select, method = "pearson")
head(cor_df)
cor_df %>% 
  as.data.frame() %>% 
  write.csv("correlation.csv")

# 读取文件
result = read.csv("correlation.csv", header = T, row.names = 1)
head(result)

head(aucs_select)
aucs_select$samples = rownames(aucs_select)
df_txndc5_data = data.frame("samples" = df$samples, "TXNDC5" = df_txndc5)
head(df_txndc5_data)

data = merge(aucs_select, df_txndc5_data, by = "samples")
head(data)
head(cor_df)

p3 = ggplot(data, aes(TXNDC5,GOBP_FIBROBLAST_PROLIFERATION)) +
  geom_point(aes(size = TXNDC5, color = TXNDC5),position = position_jitter(width = 0.2)) +
  geom_smooth(method = "lm") +
  scale_color_gradient(low = "#006fbf", high = "#ff7b00") +
  theme_bw() +
  labs(x = "TXNDC5",y = "Fibroblast proliferation")


p4 = ggplot(data, aes(TXNDC5,GOBP_FIBROBLAST_MIGRATION)) +
  geom_point(aes(size = TXNDC5, color = TXNDC5),position = position_jitter(width = 0.2)) +
  geom_smooth(method = "lm") +
  scale_color_gradient(low = "#006fbf", high = "#ff7b00") +
  theme_bw() +
  labs(x = "TXNDC5",y = "Fibroblast migration")

p5 = ggplot(data, aes(TXNDC5,GOBP_EXECUTION_PHASE_OF_APOPTOSIS)) +
  geom_point(aes(size = TXNDC5, color = TXNDC5),position = position_jitter(width = 0.2)) +
  geom_smooth(method = "lm") +
  scale_color_gradient(low = "#006fbf", high = "#ff7b00") +
  theme_bw() +
  labs(x = "TXNDC5",y = "Execution phase of apoptosis")

p6 = ggplot(data, aes(TXNDC5,GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA)) +
  geom_point(aes(size = TXNDC5, color = TXNDC5),position = position_jitter(width = 0.2)) +
  geom_smooth(method = "lm") +
  scale_color_gradient(low = "#006fbf", high = "#ff7b00") +
  theme_bw() +
  labs(x = "TXNDC5",y = "Response to transforming growth factor beta")


plot_grid(p3,p4,p5,p6, ncol = 2, labels = c("A","B","C","D"))
