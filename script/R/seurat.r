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