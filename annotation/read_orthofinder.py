def read_orthogroup_genes(file_path):
    """读取直系同源基因组文件，提取第三列的基因名称并存储到字典中"""
    
    # 读取文件，使用制表符分隔，不使用默认header
    with open(file_path, 'r') as f:
        next(f)
        lines = f.readlines()
    
    gene_dict = {}
    for line in lines:
        # 跳过空行
        if not line.strip():
            continue
            
        # 使用制表符分割行
        columns = line.strip().split('\t')
        
        # 确保至少有3列
        if len(columns) >= 3:
            # 获取第三列的基因
            genes_str = columns[2]
            if genes_str:  # 确保不是空字符串
                # 使用逗号和空格分割基因名称
                genes = genes_str.split(',')
                # 清理每个基因名称并存储到字典中
                for gene in genes:
                    gene = gene.strip()
                    if gene:  # 确保基因名称不为空
                        gene_dict[gene] = True
    
    print(f"总共找到 {len(gene_dict)} 个唯一基因")
    
    return gene_dict

def export_genes(gene_dict, output_file):
    """将基因字典导出到文件，每行一个基因"""
    with open(output_file, 'w') as f:
        for gene in sorted(gene_dict.keys()):
            f.write(f"{gene}\n")
    print(f"基因列表已保存到: {output_file}")

if __name__ == "__main__":
    input_file = r"D:\nextcloud\pd论文\data\test\Cell_type_prediction\PD-SP-cell\PD-SP-cell_Orthogroups.tsv"  # 替换为您的输入文件路径
    output_file = r"D:\nextcloud\pd论文\data\test\Cell_type_prediction\PD-SP-cell\PD-SP-cell_Orthogroups_genes.txt"  # 替换为您想要的输出文件路径
    
    # 读取基因
    gene_dict = read_orthogroup_genes(input_file)
    
    # 导出结果
    if gene_dict:
        export_genes(gene_dict, output_file)