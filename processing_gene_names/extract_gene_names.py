import pandas as pd
"""
读入marker文件,提取基因名称,marker文件格式如下
"p_val" "avg_log2FC" "pct.1" "pct.2" "p_val_adj" "cluster" "gene"
"LOC131768943" 5.98636196686601e-135 -4.71546249021902 0.016 0.352 3.21234169503997e-130 "0" "LOC131768943"
"LOC113677177" 3.38178484945174e-130 -3.41989240052076 0.02 0.354 1.8146995680643e-125 "0" "LOC113677177"
"transcript-HQ-P2-transcript13875/f192p0/1996" 2.35948385593556e-124 -2.26149823236855 0.025 0.355 1.26612263193358e-119 "0" "transcript-HQ-P2-transcript13875/f192p0/1996"
"LOC131789709" 1.9525717479567e-84 -5.2709441592786 0.009 0.225 1.04776952567104e-79 "0" "LOC131789709"
"LOC113675085" 3.78074822616779e-83 -4.3399188101396 0.01 0.225 2.0287873056439e-78 "0" "LOC113675085"
"LOC131783707" 5.03417612968189e-80 -4.85592570965736 0.007 0.214 2.7013892529486e-75 "0" "LOC131783707"

"""

marker_file = r"D:\nextcloud\pd论文\data\cluster-marker\allmarkers.tsv"
marker_df = pd.read_csv(marker_file,sep=" ")
print(marker_df.head())

# 创建列表来存储基因名称
gene_list = []
loc111_genes = []

# 遍历gene列，将基因名称添加到列表中
for gene in marker_df['gene'].unique():
    # 去除引号，如果存在的话
    gene_clean = gene.strip('"')
    gene_list.append(gene_clean)
    if gene_clean.startswith("LOC111"):
        loc111_genes.append(gene_clean)

output_file = r"D:\nextcloud\pd论文\data\cluster-marker\gene_list.tsv"
with open(output_file, 'w') as f:
    for gene in gene_list:
        f.write(f"{gene}\n")
print(f"已提取 {len(gene_list)} 个唯一基因名称并保存到 {output_file}")

output_file = r"D:\nextcloud\pd论文\data\cluster-marker\loc111_genes.tsv"
with open(output_file, 'w') as f:
    for gene in loc111_genes:
        f.write(f"{gene}\n")
print(f"已提取 {len(loc111_genes)} 个LOC111基因名称并保存到 {output_file}")