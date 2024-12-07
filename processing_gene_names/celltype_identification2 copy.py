import pandas as pd
from pathlib import Path
from itertools import product
import logging
from datetime import datetime



def read_file_with_encodings(file_path, **kwargs):
    """通用的文件读取函数，自动处理编码问题"""
    encodings = ['utf-8', 'gbk', 'gb18030']
    for encoding in encodings:
        try:
            return pd.read_csv(file_path, encoding=encoding, **kwargs)
        except UnicodeDecodeError:
            continue
    raise UnicodeDecodeError(f"无法用{encodings}中的编码读取文件：{file_path}")

def process_marker_file1(file_path, gene_to_protein, gene_prefix='LOC111', min_log2fc=0, max_pval=0.05):
    """处理第一种格式的marker文件"""
    marker_df = read_file_with_encodings(file_path, sep='\t')
    proteins = []
    
    for index, row in marker_df.iterrows():
        gene = row.iloc[0]
        avg_log2fc = row['avg_log2FC']
        p_val = row['p_val']
        
        if (gene.startswith(gene_prefix) and 
            avg_log2fc > min_log2fc and 
            p_val < max_pval and 
            gene in gene_to_protein and 
            gene_to_protein[gene] is not None):
            protein = gene_to_protein[gene].replace('.', '_')
            proteins.append(protein)
    
    return proteins

def process_marker_file2(file_path, all_marker2_files):
    """处理第二种格式的marker文件，只保留特异性表达的基因"""
    # 读取当前文件
    current_markers = read_file_with_encodings(file_path, sep=" ", quotechar='"')
    current_proteins = set()
    
    # 读取所有其他marker2文件的基因
    other_proteins = set()
    for other_file in all_marker2_files:
        if other_file != file_path:  # 跳过当前文件
            other_markers = read_file_with_encodings(other_file, sep=" ", quotechar='"')
            for gene in other_markers.index:
                if (gene.startswith('Spis-XP-') and 
                    other_markers.loc[gene, 'avg_log2FC'] > 0 and 
                    other_markers.loc[gene, 'p_val'] < 0.05):
                    protein = gene.replace('Spis-XP-', 'XP_').replace('-', '_')
                    other_proteins.add(protein)
    
    # 只保留在当前文件中特异表达的基因
    for gene in current_markers.index:
        avg_log2fc = current_markers.loc[gene, 'avg_log2FC']
        p_val = current_markers.loc[gene, 'p_val']
        
        if (gene.startswith('Spis-XP-') and 
            avg_log2fc > 0 and 
            p_val < 0.05):
            protein = gene.replace('Spis-XP-', 'XP_').replace('-', '_')
            if protein not in other_proteins:  # 只添加不在其他文件中出现的蛋白
                current_proteins.add(protein)
    
    return list(current_proteins)

def setup_logging(output_dir):
    """设置日志"""
    log_file = output_dir / f"comparison_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def save_comparison_results(comparison, output_dir):
    """保存每对比较的结果"""
    # 创建结果文件名
    result_file = output_dir / f"{comparison['marker1_file']}_vs_{comparison['marker2_file']}.tsv"
    
    # 创建结果数据框
    result_df = pd.DataFrame({
        'protein_id': sorted(list(comparison['overlapping_proteins'])),
        'source': 'overlapping'
    })
    
    # 添加只在marker1中出现的蛋白
    marker1_only = comparison['marker1_proteins'] - comparison['overlapping_proteins']
    if marker1_only:
        result_df = pd.concat([result_df, pd.DataFrame({
            'protein_id': sorted(list(marker1_only)),
            'source': 'marker1_only'
        })])
    
    # 添加只在marker2中出现的蛋白
    marker2_only = comparison['marker2_proteins'] - comparison['overlapping_proteins']
    if marker2_only:
        result_df = pd.concat([result_df, pd.DataFrame({
            'protein_id': sorted(list(marker2_only)),
            'source': 'marker2_only'
        })])
    
    # 保存结果
    result_df.to_csv(result_file, sep='\t', index=False)
    return result_file

def compare_marker_pairs(marker1_files, marker2_files, mapping_file, output_dir, 
                        gene_prefix='LOC111', min_log2fc=0, max_pval=0.05):
    """分别比较每对marker文件之间的重叠情况"""
    # 创建输出目录
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 设置日志
    logger = setup_logging(output_dir)
    
    logger.info(f"开始比较分析...")
    logger.info(f"使用参数: gene_prefix={gene_prefix}, min_log2fc={min_log2fc}, max_pval={max_pval}")
    
    # 读取mapping文件
    logger.info(f"读取mapping文件: {mapping_file}")
    mapping_df = read_file_with_encodings(mapping_file, sep='\t')
    gene_to_protein = dict(zip(mapping_df['gene_id'], mapping_df['protein_id']))
    gene_to_protein = {k: v if v != '未找到对应的protein_id' else None 
                      for k, v in gene_to_protein.items()}
    
    # 存储��个文件的处理结果
    marker1_results = {}
    marker2_results = {}
    
    # 处理所有的marker1文件
    for file in marker1_files:
        logger.info(f"处理marker1文件: {file}")
        proteins = process_marker_file1(file, gene_to_protein, 
                                      gene_prefix=gene_prefix,
                                      min_log2fc=min_log2fc,
                                      max_pval=max_pval)
        marker1_results[file] = set(proteins)
        logger.info(f"找到 {len(proteins)} 个有效蛋白ID")
    
    # 处理所有的marker2文件
    for file in marker2_files:
        logger.info(f"处理marker2文件: {file}")
        proteins = process_marker_file2(file, marker2_files)  # 传入所有marker2文件列表
        marker2_results[file] = set(proteins)
        logger.info(f"找到 {len(proteins)} 个特异性蛋白ID")
    
    # 比较每一对文件
    comparisons = []
    for marker1_file, marker2_file in product(marker1_files, marker2_files):
        marker1_name = Path(marker1_file).stem
        marker2_name = Path(marker2_file).stem
        
        logger.info(f"\n比较 {marker1_name} 和 {marker2_name}")
        
        proteins1 = marker1_results[marker1_file]
        proteins2 = marker2_results[marker2_file]
        overlapping = proteins1 & proteins2
        overlap_ratio = len(overlapping) / len(proteins1) if proteins1 else 0
        
        comparison = {
            'marker1_file': marker1_name,
            'marker2_file': marker2_name,
            'marker1_proteins': proteins1,
            'marker2_proteins': proteins2,
            'overlapping_proteins': overlapping,
            'overlap_ratio': overlap_ratio
        }
        
        # 保存这对比较的结果
        result_file = save_comparison_results(comparison, output_dir)
        logger.info(f"结果已保存到: {result_file}")
        
        # 记录统计信息
        logger.info(f"第一个文件中的有效蛋白ID数量：{len(proteins1)}")
        logger.info(f"第二个文件中的蛋白ID数量：{len(proteins2)}")
        logger.info(f"重复的蛋白ID数量：{len(overlapping)}")
        logger.info(f"重复比例：{overlap_ratio:.2%}")
        
        comparisons.append(comparison)
    
    logger.info("分析完成！")
    return comparisons

# 使用示例
if __name__ == "__main__":
    # 定义文件路径
    marker1_files = [
        r"D:\nextcloud\pd论文\data\cluster-marker\cluster0_markers.tsv",
        r"D:\nextcloud\pd论文\data\cluster-marker\cluster1_markers.tsv",
        r"D:\nextcloud\pd论文\data\cluster-marker\cluster2_markers.tsv",
        r"D:\nextcloud\pd论文\data\cluster-marker\cluster3_markers.tsv",
        r"D:\nextcloud\pd论文\data\cluster-marker\cluster6_markers.tsv",
        r"D:\nextcloud\pd论文\data\cluster-marker\cluster7_markers.tsv",
        # 添加更多第一类marker文件...
    ]
    
    marker2_files = [
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\cnidocyte-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\digestive_filaments-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\neuron-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\mitotic-host-cells-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\alga-hosting-cells-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\calicoblast-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\epidermis-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\gastrodermis-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\germline-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\gland-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\immune-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\unknown-markers.tsv"
    ]
    mapping_file = r"D:\nextcloud\pd论文\data\cluster-marker\loc111_gene_mapping.tsv"
    output_dir = r"D:\nextcloud\pd论文\result\comparison-result"

    # 运行比较
    results = compare_marker_pairs(marker1_files, marker2_files, mapping_file, output_dir)