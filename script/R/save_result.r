library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
# 添加中文字体支持
library(showtext)
library(pheatmap)
showtext_auto()  # 自动启用showtext
font_add("SimSun", "SimSun")  # 添加宋体
# 1. 数据读取和预处理
merge.data <- Read10X(data.dir = "C:/Users/伯言/Desktop/test001",gene.column = 1)
merge <- CreateSeuratObject(counts = merge.data,project = "merge")



# 1. 数据预处理和聚类（保持不变）
# 2. 数据标准化和降维
merge <- NormalizeData(merge)
all.genes <- rownames(merge)
merge <- ScaleData(merge, features = all.genes)
merge <- RunPCA(merge, features = all.genes)

# 3. UMAP降维和聚类
merge <- RunUMAP(merge, reduction = "pca", dims = 1:15)
merge <- FindNeighbors(merge, dims = 1:15)
merge <- FindClusters(merge, resolution = 0.1)

# 4. 寻找标志基因
markers <- FindAllMarkers(merge, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.25)


# 1. 获取每个聚类的细胞数量
cell_counts <- table(Idents(merge))

# 2. 获取每个聚类的基因数量
# 2. 获取每个聚类的基因数量（修改后的代码）
# 1. 初始化变量
clusters <- unique(Idents(merge))  # 获取所有聚类
cell_counts <- table(Idents(merge))  # 获取每个聚类的细胞数量

# 2. 获取每个聚类的基因数量
gene_counts <- numeric(length(clusters))
expression_matrix <- GetAssayData(merge, slot = "counts")  # 使用标准化后的数据

for(i in seq_along(clusters)) {
  # 获取当前聚类的细胞
  cluster_cells <- WhichCells(merge, idents = clusters[i])
  # 获取当前聚类的表达矩阵
  cluster_data <- expression_matrix[, cluster_cells]
  # 计算在该聚类中至少有一个细胞表达的基因数量
  gene_counts[i] <- sum(rowSums(cluster_data > 0) > 0)
}

# 3. 获取每个聚类的标记基因数量
marker_counts <- table(markers$cluster)

# 4. 创建结果数据框
cluster_stats <- data.frame(
  Cluster = clusters,
  Cell_Count = as.numeric(cell_counts),
  Gene_Count = gene_counts,
  Marker_Gene_Count = as.numeric(marker_counts)
)

# 5. 添加总计行
total_row <- data.frame(
  Cluster = "Total",
  Cell_Count = sum(cluster_stats$Cell_Count),
  Gene_Count = sum(cluster_stats$Gene_Count),
  Marker_Gene_Count = sum(cluster_stats$Marker_Gene_Count)
)

cluster_stats <- rbind(cluster_stats, total_row)

# 6. 保存为tsv文件
write.table(cluster_stats, 
            file = "cluster_statistics4.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# 获取每个聚类的标记基因统计
cluster_gene_stats <- markers %>%
  group_by(cluster) %>%
  summarise(
    total_genes = n(),  # 所有差异基因
    up_genes = sum(avg_log2FC > 0),  # 上调基因
    down_genes = sum(avg_log2FC < 0)  # 下调基因
  )
# 保存详细统计结果
write.table(cluster_gene_stats,
            file = "cluster_gene_statistics2.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
