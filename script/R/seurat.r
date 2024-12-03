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

## qc

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
markers <- FindAllMarkers(merge, only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.25)


# 4.1 为每个聚类保存标志基因
clusters <- unique(markers$cluster)
# 创建输出目录
dir.create("out/markers_by_cluster",recursive = TRUE,showWarnings = FALSE)
# 按照聚类分别保存
for(cluster in clusters){
  cluster_markers <- markers %>%
    filter(cluster == !!cluster) %>%
    arrange(desc(avg_log2FC)) #按照差异倍数来排序
  # 构建输出名
  output_dir <- sprintf("out/markers_by_cluster/cluster%s_markers1.tsv",cluster)
  # 保存文件
  write.table(cluster_markers,
              file = output_dir,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  # 打印进度信息
  print(sprintf("已经保存cluster %s 的marker基因，共%d个基因",cluster,nrow(cluster_markers)))
  
}

write.table(markers, 
            file = "out/markers_by_cluster/all_markers.tsv",  # 输出文件名
            sep = "\t",                # 使用tab作为分隔符
            quote = FALSE,             # 不使用引号
            row.names = FALSE)         # 不包含行名

# 5. 基本统计
# 查看聚类数量和细胞数
print("每个聚类的细胞数：")
print(table(merge$seurat_clusters))

# 查看每个聚类的标志基因数量
print("每个聚类的标志基因数：")
print(markers %>%
        group_by(cluster) %>%
        summarise(gene_count = n()) %>%
        arrange(cluster))


# 6. 核心可视化
# UMAP聚类图
umap_plot <- DimPlot(merge, reduction = "umap", label = TRUE, repel = TRUE)
ggsave("umap_clusters.pdf", umap_plot, width = 10, height = 8)

# 1. 先获取基础数据
cluster_cells <- table(merge$seurat_clusters)
cluster_genes <- markers %>%
  group_by(cluster) %>%
  summarise(gene_count = n()) %>%
  mutate(cluster = as.numeric(as.character(cluster)))  # 转换cluster为数值类型

# 2. 创建数据框
cluster_stats <- data.frame(
  cluster = as.numeric(names(cluster_cells)),
  cell_count = as.numeric(cluster_cells)
) %>%
  left_join(cluster_genes, by = "cluster")

# 3. 绘制统计图
# 细胞数量分布
p1 <- ggplot(cluster_stats, aes(x = factor(cluster), y = cell_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "每个聚类的细胞数量",  # 显示总细胞数
       x = "聚类编号",
       y = "细胞数量") +
  theme(plot.title = element_text(hjust = 0.5))

# 基因数量分布
p2 <- ggplot(cluster_stats, aes(x = factor(cluster), y = gene_count)) +
  geom_bar(stat = "identity", fill = "coral") +
  theme_minimal() +
  labs(title = "每个聚类的基因数量", # 显示总基因数
       x = "聚类编号",
       y = "基因数量") +
  theme(plot.title = element_text(hjust = 0.5))

# 组合图形
combined_plot <- p1 / p2

# 保存统计图
ggsave("cluster_statistics.pdf", 
       combined_plot, 
       width = 6, 
       height = 6)


# 1. 准备数据时进行排序
# 按聚类对细胞进行排序
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%  # 每个cluster选择top 10的基因
  pull(gene)  # 提取基因名

# 2. 准备热图数据
expr_matrix <- GetAssayData(merge, slot = "scale.data")[top_markers, cell_order]

# 2. 创建排序后的注释数据
annotation_col <- data.frame(
  Cluster = merge$seurat_clusters[cell_order],
  row.names = colnames(merge)[cell_order]
)
ann_colors <- list(
  Cluster = c(
    "0" = "#E41A1C",  # 红色
    "1" = "#377EB8",  # 蓝色
    "2" = "#4DAF4A",  # 绿色
    "3" = "#984EA3",  # 紫色
    "4" = "#FF7F00",  # 橙色
    "5" = "#FFFF33",  # 黄色
    "6" = "#A65628",  # 棕色
    "7" = "#F781BF",  # 粉色
    "8" = "#999999"   # 灰色
    
  )
)

# 3. 在pheatmap中禁用列聚类
pdf("marker_genes_heatmap1.pdf", width = 12, height = 8)
pheatmap(
  expr_matrix,
  annotation_col = annotation_col,
  ann_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_cols = FALSE,  # 关键：禁用列聚类
  cluster_rows = FALSE,
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#800080", "#000000", "#FFD700"))(100),
  gaps_col = cumsum(table(merge$seurat_clusters[cell_order])),  # 添加分隔线
  treeheight_row = 20,
  treeheight_col = 0,  # 由于不聚类，可以移除列树状图
  border_color = NA,
  fontsize = 8,
  fontsize_row = 6,
  legend = TRUE,
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2", "-1", "0", "1", "2"),
  main = "Marker Genes Expression Heatmap",
  fontfamily = "SimSun"
)
dev.off()


plants_genes <- grepl("SPIL2461", rownames(merge))
sum(plants_genes)

all.genes[plants_genes]

all_plants_genes <- all.genes[plants_genes]

# 创建一个新的"合成基因"，将所有植物基因的表达量相加
plant_sum <- colSums(GetAssayData(merge)[all_plants_genes,])

# 将这个新的"合成基因"添加到meta.data中
merge$plant_sum <- plant_sum

# 绘制这个合成基因的表达图
p <- FeaturePlot(merge, 
                 features = "plant_sum",
                 reduction = "umap")
# 保存图像
ggsave("plant_sum_expression.pdf", 
       plot = p, 
       width = 10, 
       height = 8, 
       dpi = 300) # 分辨率

plants_genes <- grepl("SPIL2461", rownames(merge))
all.genes[plants_genes]



cluster.markers <- FindMarkers(merge,ident.1=c(4, 5,8,9))


write.table(cluster.markers,
            file = "combined_clusters_vs_459_markers2.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

####      聚类中的基因数量与marker基因数量对照



# 统计每个聚类的标记基因与总基因的对照表
# 获取每个聚类的总基因数
total_genes_per_cluster <- sapply(levels(Idents(merge)), function(cluster) {
  cluster_cells <- WhichCells(merge, idents = cluster)
  cluster_data <- GetAssayData(merge, slot = "counts")[, cluster_cells]
  sum(rowSums(cluster_data) > 0)
})

# 获取每个聚类的标记基因数
marker_genes_per_cluster <- markers %>%
  group_by(cluster) %>%
  summarise(marker_count = n())

# 创建对照表
comparison_table <- data.frame(
  Cluster = levels(Idents(merge)),
  Total_Genes = total_genes_per_cluster,
  Marker_Genes = marker_genes_per_cluster$marker_count,
  Percentage = round(marker_genes_per_cluster$marker_count / total_genes_per_cluster * 100, 2)
)

# 添加总计行
total_row <- data.frame(
  Cluster = "Total",
  Total_Genes = sum(total_genes_per_cluster),
  Marker_Genes = sum(marker_genes_per_cluster$marker_count),
  Percentage = round(sum(marker_genes_per_cluster$marker_count) / sum(total_genes_per_cluster) * 100, 2)
)

comparison_table <- rbind(comparison_table, total_row)

# 保存对照表
write.table(comparison_table,
            file = "marker_vs_total_genes_comparison.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)



# 获取共生细胞聚类的差异基因详情
symbiont_clusters <- c(4, 5, 8, 9)
symbiont_markers <- cluster.markers

# 统计共生细胞数量
symbiont_cells <- sum(merge$seurat_clusters %in% symbiont_clusters)
# 假设 symbiont_markers 是你想要保存的数据框
write.table(symbiont_markers, file="symbiont_markers.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)




# 寻找4589合并后的标记基因
cluster2.marker <- FindMarkers(merge,ident.1 = c(4,5,8,9))
head(cluster2.marker,n=10)