import pandas as pd
import glob

# 获取所有的marker文件
marker_files = glob.glob("data/4589/cluster*_markers.tsv")

# 处理每个文件
for file in sorted(marker_files):
    # 读取TSV文件
    df = pd.read_csv(file, sep='\t')
    
    # 只保留gene列
    genes = df['gene'].tolist()
    
    # 生成输出文件名
    output_file = file.replace('_markers.tsv', '_genes.txt')
    
    # 保存到新文件
    with open(output_file, 'w') as f:
        for gene in genes:
            f.write(f"{gene}\n")
