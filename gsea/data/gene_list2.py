import pandas as pd

# 读取数据
data = pd.read_csv(r'D:\nextcloud\pd论文\data\cluster-marker\cluster0_markers.tsv',sep='\t')  # 假设文件名为gene_data.csv
data.rename(columns={'Unnamed: 0': 'gene_name'},inplace=True)


print(data.columns)  # 添加此行以查看数据框的列名

# 根据log2FC的正负值分组
positive_genes = data[data['avg_log2FC'] > 0]['gene_name'].str.replace('-','_')  # 正值基因
negative_genes = data[data['avg_log2FC'] < 0]['gene_name'].str.replace('-','_')  # 负值基因

# 保存到文件
positive_genes.to_csv(r'D:\nextcloud\pd论文\data\genelist\Cluster0\positive_genes.txt', index=False, header=False)
negative_genes.to_csv(r'D:\nextcloud\pd论文\data\genelist\Cluster0\negative_genes.txt', index=False, header=False)


