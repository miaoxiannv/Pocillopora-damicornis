import pandas as pd

"""
读取文件格式为，只存储protein_id，不存储gene_id
gene_id	protein_id
LOC111325671	XP_022785263.1
LOC111326254	XP_022785952.1
LOC111326848	XP_022786664.1
LOC111339103	未找到对应的protein_id
LOC111343026	XP_022805896.1
LOC111322693	XP_022781572.1
LOC111320871	XP_022779335.1
LOC111340891	XP_022803549.1
LOC111327082	XP_022786918.1
LOC111320821	XP_022779255.1
LOC111328748	XP_022788977.1
"""

# 读取映射文件，尝试不同的编码方式
try:
    # 首先尝试 utf-8
    mapping_file = r"D:\nextcloud\pd论文\data\cluster-marker\loc111_gene_mapping.tsv"
    gene_protein_map = pd.read_csv(mapping_file, sep='\t')
except UnicodeDecodeError:
    try:
        # 如果 utf-8 失败，尝试 GBK
        gene_protein_map = pd.read_csv(mapping_file, sep='\t', encoding='gbk')
    except UnicodeDecodeError:
        # 如果 GBK 也失败，尝试 gb18030
        gene_protein_map = pd.read_csv(mapping_file, sep='\t', encoding='gb18030')

# 打印基本信息
total_genes = len(gene_protein_map)
found_proteins = len(gene_protein_map[gene_protein_map['protein_id'] != '��找到对应的protein_id'])
print(f"总共读取了 {total_genes} 个基因")
print(f"其中有 {found_proteins} 个基因找到了对应的protein_id")
print(f"有 {total_genes - found_proteins} 个基因未找到对应的protein_id")
"""
读取marker基因文件，文件格式为
"p_val" "avg_log2FC" "pct.1" "pct.2" "p_val_adj" "cluster" "gene"
"Spis-XP-022798554-1" 0 9.74859655423309 0.878 0.015 0 "" "Spis-XP-022798554-1"
"Spis-XP-022781821-1" 0 9.50037024885467 0.732 0.008 0 "" "Spis-XP-022781821-1"
"Spis-XP-022808637-1" 0 9.0280294531793 0.439 0.003 0 "" "Spis-XP-022808637-1"
"orphan-peak-3033" 0 8.13707132316073 0.415 0.005 0 "" "orphan-peak-3033"
"Spis-XP-022798199-1" 0 7.99327225918787 0.39 0.005 0 "" "Spis-XP-022798199-1"
"Spis-XP-022780446-1" 0 7.53360173604392 0.293 0.003 0 "" "Spis-XP-022780446-1"
"Spis-XP-022781196-1" 2.78359796592752e-307 7.96039948900419 0.415 0.008 1.04050891966371e-302 "" "Spis-XP-022781196-1"
"Spis-XP-022804742-1" 4.73834683563687e-297 8.87662734596142 0.28 0.003 1.77119404716106e-292 "" "Spis-XP-022804742-1"
"Spis-LOC111338138" 7.31710161075513e-276 7.48449944957496 0.244 0.003 2.73513258210027e-271 "" "Spis-LOC111338138"
"Spis-XP-022809999-1" 2.99599701211422e-261 7.67891932175126 0.256 0.003 1.11990368312829e-256 "" "Spis-XP-022809999-1"
"""

# 读取marker基因文件，同样处理编码问题
try:
    marker_file = r"D:\nextcloud\pd论文\data\Cell-cluster-marker\allmarkers.tsv"
    marker_df = pd.read_csv(marker_file, sep=" ")
except UnicodeDecodeError:
    try:
        marker_df = pd.read_csv(marker_file, sep=" ", encoding='gbk')
    except UnicodeDecodeError:
        marker_df = pd.read_csv(marker_file, sep=" ", encoding='gb18030')

# 提取基因名称并清理格式
marker_genes = marker_df['gene'].str.strip('"').unique()

# 处理基因名称，去除"Spis-"前缀
cleaned_genes = []
for gene in marker_genes:
    if gene.startswith('Spis-'):
        gene = gene[5:]  # 去除"Spis-"前缀
        gene = gene.replace('-', '_')
    cleaned_genes.append(gene)

# 统计基因类型
xp_genes = sum(1 for gene in cleaned_genes if gene.startswith('XP_'))
loc_genes = sum(1 for gene in cleaned_genes if gene.startswith('LOC'))
other_genes = len(cleaned_genes) - xp_genes - loc_genes

print(f"总共读取了 {len(cleaned_genes)} 个marker基因")
print(f"其中：")
print(f"XP开头的基因数量：{xp_genes}")
print(f"LOC开头的基因数量：{loc_genes}")
print(f"其他类型的基因数量：{other_genes}")

# 获取有效的protein_id列表（排除未找到的）
valid_proteins = gene_protein_map[gene_protein_map['protein_id'] != '未找到对应的protein_id']['protein_id'].tolist()
valid_proteins = [p.replace('.', '_') for p in valid_proteins]
# 将XP基因提取出来（之前已经去除了Spis-前缀）
marker_xp_genes = [gene for gene in cleaned_genes if gene.startswith('XP_')]

# 查找重复的基因
overlapping_genes = set(valid_proteins) & set(marker_xp_genes)

# 打印重叠信息
print("\n基因重叠分析：")
print(f"蛋白质ID数量：{len(valid_proteins)}")
print(f"Marker中XP基因数量：{len(marker_xp_genes)}")
print(f"重叠的基因数量：{len(overlapping_genes)}")

# 如果需要，可以保存重叠的基因列表
#overlap_output = "overlapping_genes.tsv"
# with open(overlap_output, 'w') as f:
#   for gene in sorted(overlapping_genes):
#        f.write(f"{gene}\n")

#print(f"重叠的基因列表已保存到：{overlap_output}")

# 打印前几个重叠的基因作为示例
print("\n重叠基因示例（前5个）：")
for gene in sorted(list(overlapping_genes))[:5]:
    print(gene)
