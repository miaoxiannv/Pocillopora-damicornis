import os
from pathlib import Path

import pandas as pd

"""
读取marker文件将正负marker分开
marker文件格式：
	p_val	avg_log2FC	pct.1	pct.2	p_val_adj
LOC131768943	5.98636196686601e-135	-4.71546249021902	0.016	0.352	3.21234169503997e-130
LOC113677177	3.38178484945174e-130	-3.41989240052076	0.02	0.354	1.8146995680643e-125
transcript-HQ-P2-transcript13875/f192p0/1996	2.35948385593556e-124	-2.26149823236855	0.025	0.355	1.26612263193358e-119
LOC131789709	1.9525717479567e-84	-5.2709441592786	0.009	0.225	1.04776952567104e-79
LOC113675085	3.78074822616779e-83	-4.3399188101396	0.01	0.225	2.0287873056439e-78
LOC131783707	5.03417612968189e-80	-4.85592570965736	0.007	0.214	2.7013892529486e-75
transcript-HQ-P2-transcript20/f2p0/8213	1.07825783316943e-75	-3.36122210049257	0.01	0.211	5.78603935857048e-71
LOC113673226	2.84635829214506e-72	-2.98720465527355	0.012	0.209	1.52738432314796e-67
transcript-HQ-OA1-1SP-transcript21/f3p0/8153	4.94409833349172e-71	-4.15686571999901	0.008	0.195	2.65305260673499e-66
LOC131769462	2.50178660230857e-70	-5.0675436273289	0.022	0.224	1.3424837086648e-65
LOC113679383	2.62554746056298e-70	-4.14358064801303	0.025	0.23	1.4088950228127e-65
"""


def process_multiple_markers(input_dir, output_dir):
    """
    批量处理多个marker文件
    
    Parameters:
    input_dir (str): 包含所有marker文件的输入目录
    output_dir (str): 输出文件的目录
    
    Returns:
    dict: 包含每个文件处理统计信息的字典
    """
    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)

    # 存储所有文件的统计信息
    all_stats = {}

    # 处理所有的.tsv文件
    for file in Path(input_dir).glob("*markers.tsv"):
        cluster_name = file.stem  # 获取文件名（不含扩展名）

        # 设置输出文件路径
        output_positive = os.path.join(output_dir, f"{cluster_name}_positive.txt")
        output_negative = os.path.join(output_dir, f"{cluster_name}_negative.txt")

        # 读取并处理文件
        try:
            df = pd.read_csv(file, sep='\t')

            # 首先过滤 p_val < 0.05 的数据
            df = df[df['p_val'] < 0.05]
            if 'Unnamed: 0' in df.columns:
                df.set_index(df['Unnamed: 0'], inplace=True)
                df = df.drop('Unnamed: 0', axis=1)
            # 分离正负marker
            positive_markers = df[df['avg_log2FC'] > 0]
            negative_markers = df[df['avg_log2FC'] < 0]
            # 过滤掉表头并转换为字符串列表
            positive_genes = [str(x) for x in positive_markers.index if str(x) != 'Unnamed: 0']
            negative_genes = [str(x) for x in negative_markers.index if str(x) != 'Unnamed: 0']
            # 排序
            positive_markers = positive_markers.sort_values('avg_log2FC', ascending=False)
            negative_markers = negative_markers.sort_values('avg_log2FC', ascending=True)
            # 保存结果
            #positive_markers.to_csv(output_positive, sep='\t', index=True)
            #negative_markers.to_csv(output_negative, sep='\t', index=True)
            # 保存基因名称，确保没有表头
            with open(output_positive, 'w', encoding='utf-8') as f:
                f.write('\n'.join(positive_genes))
            
            with open(output_negative, 'w', encoding='utf-8') as f:
                f.write('\n'.join(negative_genes))
            # 记录统计信息
            all_stats[cluster_name] = {
                'total_markers': len(df),
                'positive_markers': len(positive_markers),
                'negative_markers': len(negative_markers)
            }

            print(f"Processed {cluster_name}: {len(df)} total markers")

        except Exception as e:
            print(f"Error processing {file}: {str(e)}")
            all_stats[cluster_name] = {'error': str(e)}

    return all_stats


# 使用示例
input_dir = r"D:\nextcloud\pd论文\data\cluster-marker"
output_dir = r"D:\nextcloud\pd论文\data\Cluster-pos-nag-genlist"
stats = process_multiple_markers(input_dir, output_dir)

# 打印统计信息
for cluster, stat in stats.items():
    if 'error' not in stat:
        print(f"\nCluster: {cluster}")
        print(f"Total markers: {stat['total_markers']}")
        print(f"Positive markers: {stat['positive_markers']}")
        print(f"Negative markers: {stat['negative_markers']}")
