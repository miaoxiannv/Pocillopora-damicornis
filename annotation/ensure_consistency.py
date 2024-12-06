def read_genes_from_file(file_path):
    """从文件中读取基因名称，返回一个集合"""
    with open(file_path, 'r') as f:
        genes = {line.strip() for line in f if line.strip()}
    return genes

def find_overlapping_genes(file1, file2):
    """查找两个文件中重叠的基因"""
    genes1 = read_genes_from_file(file1)
    genes2 = read_genes_from_file(file2)
    
    # 找到重叠的基因
    overlapping_genes = genes1.intersection(genes2)
    
    print(f"文件1中基因数: {len(genes1)}")
    print(f"文件2中基因数: {len(genes2)}")
    print(f"重叠的基因数: {len(overlapping_genes)}")
    
    return overlapping_genes

if __name__ == "__main__":
    file1 = r"D:\nextcloud\pd论文\data\test\Cell_type_prediction\PD-SP-cell\PD-SP-cell_Orthogroups_genes.txt"  # 替换为第一个文件的路径
    file2 = r"D:\nextcloud\pd论文\data\test\Cell_type_prediction\PD-SP-cell\feature_names.tsv"  # 替换为第二个文件的路径
    
    # 查找重叠的基因
    overlapping_genes = find_overlapping_genes(file1, file2)
    
    # 输出重叠的基因
    if overlapping_genes:
        print("重叠的基因数量:")
        print(len(overlapping_genes))
