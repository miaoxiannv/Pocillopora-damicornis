import pandas as pd

"""
读取mapping文件格式如下
gene_id	protein_id
LOC111326392	未找到对应的protein_id
LOC111326177	XP_022785962.1
...
"""

# 读取第一个mapping文件
try:
    mapping_file = r"D:\nextcloud\pd论文\data\cluster-marker\ncbi_gene_mapping.tsv"  # 请替换为实际的文件路径
    mapping_df = pd.read_csv(mapping_file, sep='\t')
except UnicodeDecodeError:
    try:
        mapping_df = pd.read_csv(mapping_file, sep='\t', encoding='gbk')
    except UnicodeDecodeError:
        mapping_df = pd.read_csv(mapping_file, sep='\t', encoding='gb18030')

# 将第一个文件的protein_id中的点号替换为下划线
mapping_df['protein_id'] = mapping_df['protein_id'].apply(lambda x: x.replace('.', '_') if isinstance(x, str) else x)

# 获取第一个文件的有效蛋白ID列表（排除未找到的情况）
proteins1 = [p for p in mapping_df['protein_id'] if p != '未找到对应的protein_id']

"""
读取文件格式如下
alga-hosting_cells calicoblast cnidocyte ...
Spis10006_1 1.00647288442868 1.02293552443941 ...
Spis10012_1 1.00018842005584 1.02960418595757 ...
...
"""

# 读取文件
try:
    file_path = r"D:\nextcloud\pd论文\论文素材\cell_sp论文\g2mysjfp52-3\Spis_adult_cell_type_gene_FC.tsv\Spis_coral_cell_type_gene_FC.tsv"  # 请替换为实际的文件路径
    # 只读取第一列，将第一列设为索引
    genes_df = pd.read_csv(file_path, sep='\t')
except UnicodeDecodeError:
    try:
        genes_df = pd.read_csv(file_path, sep='\t', usecols=[0], encoding='gbk')
    except UnicodeDecodeError:
        genes_df = pd.read_csv(file_path, sep='\t', usecols=[0], encoding='gb18030')
# 处理第二个文件的基因名（去除Spis-前缀）

proteins2 = genes_df.index.tolist()
cleaned_proteins = [protein.replace('Spis_', '') if protein.startswith('Spis_') else protein for protein in proteins2]
# 找出重复的蛋白ID
overlapping_proteins = set(proteins1) & set(cleaned_proteins)

# 打印统计信息
print("\n蛋白ID统计：")
print(f"NCBI物种的有效蛋白ID数量：{len(proteins1)}")
print(f"文章中的蛋白ID数量：{len(proteins2)}")
print(f"重复的蛋白ID数量：{len(overlapping_proteins)}")
print(f"重复比例（相对于NCBI物种）：{len(overlapping_proteins)/len(proteins1):.2%}")

# 打印重复蛋白ID的示例
print("\n重复蛋白ID示例（前5个）：")
for protein in sorted(list(overlapping_proteins))[:5]:
    print(protein)