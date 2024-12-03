import os

from cell_type_processor import CellTypeProcessor
import csv
import pandas as pd

def main(DATA_PAIRS,n):
    # 初始化处理器
    processor = CellTypeProcessor()

    # 处理所有物种对
    results = processor.process_species_pairs(DATA_PAIRS)

    # 收集预测数据
    prediction_data = collect_prediction_data(results)


    # 转换为 DataFrame
    df = pd.DataFrame(prediction_data, columns=['cluster', 'predicted_type', 'species_pair', 'confidence', 'features'])

    # 统计每种细胞类型的支持数量
    support_count = df.groupby(['cluster', 'predicted_type']).size().reset_index(name='support_count')

    # 计算加权置信度
    weighted_confidence = df.groupby(['cluster', 'predicted_type']).agg({'confidence': 'mean'}).reset_index()
    weighted_confidence = weighted_confidence.merge(support_count, on=['cluster', 'predicted_type'])

    # 计算加权置信度（可以根据需要调整权重）
    weighted_confidence['weighted_confidence'] = weighted_confidence['confidence'] * weighted_confidence[
        'support_count']

    # 找到每个聚类中加权置信度最高的细胞类型
    result = weighted_confidence.loc[weighted_confidence.groupby('cluster')['weighted_confidence'].idxmax()]
    top_n = weighted_confidence.groupby('cluster').apply(lambda x: x.nlargest(n, 'weighted_confidence')).reset_index(
        drop=True)

    # 输出结果
    print(result[['cluster', 'predicted_type', 'support_count', 'weighted_confidence']])
    print(top_n[['cluster', 'predicted_type', 'support_count', 'weighted_confidence']])
    top_n[['cluster', 'predicted_type', 'support_count', 'weighted_confidence']].to_csv('top_n_cells.tsv', sep='\t', index=False)



def collect_prediction_data(results):
    """收集预测数据的辅助函数

    Args:
        results (dict): 映射结果字典，格式为 {species_pair: {cell_type: [mappings]}}

    Returns:
        list: 预测数据列表，每个元素包含 cluster, predicted_type, confidence, features 等信息
    """
    predictions = []

    for species_pair, cell_type_mappings in results.items():
        for cell_type, mappings in cell_type_mappings.items():
            for mapping in mappings:  # 遍历每个映射
                # 获取最佳映射结果
                confidence = mapping['confidence']
                overlap_count = mapping['overlap_count']
                overlap_genes = mapping['overlap_genes']

                # 提取特征
                features = [
                    overlap_count,
                    confidence,
                    len(overlap_genes)
                ]

                # 创建预测数据
                prediction = {
                    'cluster': cell_type,
                    'predicted_type': mapping['sp1_cell_type'],  # 使用 sp1_cell_type
                    'species_pair': species_pair,
                    'confidence': confidence,
                    'features': features
                }
                predictions.append(prediction)
    return predictions


def save_results_to_tsv(results, output_dir):
    """将结果保存为 TSV 文件
    
    Args:
        results (dict): 映射结果字典，格式为 {species_pair: {cell_type: [mappings]}}
        output_dir (str): 输出目录
    """
    tsv_file_path = os.path.join(output_dir, 'results.tsv')

    with open(tsv_file_path, mode='w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerow(['cluster','species_pair', 'cell_type', 'confidence', 'overlap_count', 'overlap_genes'])  # 写入表头

        for cluster, species_pair, cell_type_mappings in results.items():
            for cell_type, mappings in cell_type_mappings.items():
                if mappings:  # 确保有映射结果
                    best_mapping = max(mappings, key=lambda x: x['confidence'])
                    writer.writerow([
                        species_pair,
                        cell_type,
                        best_mapping['confidence'],
                        best_mapping['overlap_count'],
                        ','.join(best_mapping['overlap_genes'])  # 将重叠基因连接为字符串
                    ])


if __name__ == "__main__":
    DATA_PAIRS = {
        'FL-SP': {
            'ortho': r'D:\nextcloud\pd论文\data\test\Cell_type_prediction\FL-SP-cell\FL-SP-cell_Orthogroups.tsv',
            'marker1': r'D:\nextcloud\pd论文\data\test\Cell_type_prediction\FL-SP-cell\SP-cell_marker.tsv',
            'marker2': r'D:\nextcloud\pd论文\data\test\Cell_type_prediction\FL-SP-cell\all_markers.tsv'
        }, 'PD-SP': {
            'ortho': r'D:\nextcloud\pd论文\data\test\Cell_type_prediction\PD-SP-cell\PD-SP-cell_Orthogroups.tsv',
            'marker1': r'D:\nextcloud\pd论文\data\test\Cell_type_prediction\PD-SP-cell\SP-cell_marker.tsv',
            'marker2': r'D:\nextcloud\pd论文\data\test\Cell_type_prediction\PD-SP-cell\all_markers.tsv'
        }
        , 'PV-SP': {
            'ortho': r'D:\nextcloud\pd论文\data\test\Cell_type_prediction\PV-SP-cell\PV-SP-cell_Orthogroups.tsv',
            'marker1': r'D:\nextcloud\pd论文\data\test\Cell_type_prediction\PV-SP-cell\SP-cell_marker.tsv',
            'marker2': r'D:\nextcloud\pd论文\data\test\Cell_type_prediction\PV-SP-cell\all_markers.tsv'
        }
        , 'SP-SP': {
            'ortho': r'D:\nextcloud\pd论文\data\test\Cell_type_prediction\SP-SP-cell\SP-SP-cell_Orthogroups.tsv',
            'marker1': r'D:\nextcloud\pd论文\data\test\Cell_type_prediction\SP-SP-cell\SP-cell_marker.tsv',
            'marker2': r'D:\nextcloud\pd论文\data\test\Cell_type_prediction\SP-SP-cell\all_markers.tsv'
        }
    }
    #输出前n个置信度较高的预测结果
    n = 3
    main(DATA_PAIRS,n)
