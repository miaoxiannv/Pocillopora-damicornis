import pandas as pd

"""
读取marker文件格式如下
p_val	avg_log2FC	pct.1	pct.2	p_val_adj
LOC131768943	5.98636196686601e-135	-4.71546249021902	0.016	0.352	3.21234169503997e-130
LOC113677177	3.38178484945174e-130	-3.41989240052076	0.02	0.354	1.8146995680643e-125
...
"""

# 读取第一个marker文件（cluster0_markers.tsv）
try:
    marker_file1 = r"D:\nextcloud\pd论文\data\cluster-marker\cluster0_markers.tsv"
    marker_df1 = pd.read_csv(marker_file1, sep='\t')
except UnicodeDecodeError:
    try:
        marker_df1 = pd.read_csv(marker_file1, sep='\t', encoding='gbk')
    except UnicodeDecodeError:
        marker_df1 = pd.read_csv(marker_file1, sep='\t', encoding='gb18030')

"""
读取mapping文件格式如下
gene_id	protein_id
LOC111325671	XP_022785263.1
LOC111326254	XP_022785952.1
...
"""
# 读取mapping文件
try:
    mapping_file = r"D:\nextcloud\pd论文\data\cluster-marker\loc111_gene_mapping.tsv"
    mapping_df = pd.read_csv(mapping_file, sep='\t')
except UnicodeDecodeError:
    try:
        mapping_df = pd.read_csv(mapping_file, sep='\t', encoding='gbk')
    except UnicodeDecodeError:
        mapping_df = pd.read_csv(mapping_file, sep='\t', encoding='gb18030')

# 创建gene_id到protein_id的映射字典，排除未找到的情况
gene_to_protein = dict(zip(
    mapping_df['gene_id'], 
    mapping_df['protein_id']
))

# 将未找到对应protein_id的基因对应的值设为None
gene_to_protein = {k: v if v != '未找到对应的protein_id' else None for k, v in gene_to_protein.items()}

"""
读取第二个marker文件，格式如下
"p_val" "avg_log2FC" "pct.1" "pct.2" "p_val_adj"
"Spis-XP-022783415-1" 0 6.1005895284779 0.839 0.042 0
"Spis-LOC111326826" 0 3.28421454913096 0.782 0.174 0
...
"""
# 读取第二个marker文件（带有Spis-前缀的文件）
try:
    marker_file2 = r"D:\nextcloud\pd论文\data\Cell-cluster-marker\calicoblast-markers.tsv"
    marker_df2 = pd.read_csv(marker_file2, sep=" ", quotechar='"')
except UnicodeDecodeError:
    try:
        marker_df2 = pd.read_csv(marker_file2, sep=" ", quotechar='"', encoding='gbk')
    except UnicodeDecodeError:
        marker_df2 = pd.read_csv(marker_file2, sep=" ", quotechar='"', encoding='gb18030')

# 从第二个marker文件中提取蛋白ID

marker2_proteins = []
for gene in marker_df2.index:
    if gene.startswith('Spis-XP-'):
        # 去除Spis-前缀，将-替换为_
        protein = gene.replace('Spis-XP-', 'XP_').replace('-', '_')
        marker2_proteins.append(protein)




marker1_proteins = []
for gene in marker_df1.iloc[:, 0]:
    if gene.startswith('LOC111') and gene in gene_to_protein and gene_to_protein[gene] is not None:
        protein = gene_to_protein[gene].replace('.', '_')
        marker1_proteins.append(protein)

# 找出重复的蛋白ID
overlapping_proteins = set(marker1_proteins) & set(marker2_proteins)


# 计算重复比例
overlap_ratio = len(overlapping_proteins) / len(marker1_proteins) if len(marker1_proteins) > 0 else 0


# 打印统计信息
print("\n蛋白ID比对统计：")
print(f"第一个marker文件中的有效蛋白ID数量：{len(marker1_proteins)}")
print(f"第二个marker文件中的蛋白ID数量：{len(marker2_proteins)}")
print(f"重复的蛋白ID数量：{len(overlapping_proteins)}")
print(f"重复比例（重复数/marker1数量）：{overlap_ratio:.2%}")


# 打印重复蛋白ID的示例
print("\n重复蛋白ID示例（前5个）：")
for protein in sorted(list(overlapping_proteins))[:5]:
    print(protein)