import os
import pandas as pd

class AnnotationCellType:
    def __init__(self, marker_file, annotation_file):
        self.marker_file = marker_file
        self.annotation_file = annotation_file
        self.marker_df = None
        self.annotation_df = None
        self.merged_df = None
    
    def read_marker_file(self):
        # 读取marker文件，包含p_val等统计信息
        self.marker_df = pd.read_csv(self.marker_file, sep='\t', header=0)
        
        # 确保必需的统计列存在
        required_columns = ['Gene_id', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj']
        for col in required_columns:
            if col not in self.marker_df.columns:
                print(f"警告: 缺少列 {col}")
                self.marker_df[col] = ''
        
        # 标准化基因ID
        self.marker_df['Gene_id'] = self.marker_df['Gene_id'].str.strip().str.lower()
        return self.marker_df
    
    def read_annotation_file(self):
        # 读取注释文件，不使用默认header
        self.annotation_df = pd.read_csv(self.annotation_file, sep='\t', header=None)
        
        # 根据文件内容设置列名
        # 第一列为Gene_id，第二列为Domain，其余列合并为Description
        columns = ['Gene_id', 'Domain']
        for i in range(2, len(self.annotation_df.columns)):
            columns.append(f'Description{i-1}')
        
        self.annotation_df.columns = columns
        
        # 标准化基因ID
        self.annotation_df['Gene_id'] = self.annotation_df['Gene_id'].str.strip().str.lower()
        
        # 清理Domain列（去除末尾的斜杠）
        self.annotation_df['Domain'] = self.annotation_df['Domain'].str.rstrip('/')
        
        return self.annotation_df
    
    def merge_data(self):
        # 合并marker信息和注释信息
        self.merged_df = pd.merge(self.marker_df, 
                                self.annotation_df,
                                on='Gene_id',
                                how='left')
        
        # 统计合并结果
        total_genes = len(self.marker_df)
        annotated_genes = self.merged_df['Domain'].notna().sum()
        print(f"总基因数: {total_genes}")
        print(f"获得注释的基因数: {annotated_genes}")
        print(f"未注释的基因数: {total_genes - annotated_genes}")
        
        return self.merged_df
    
    def export(self, out_file):
        # 确保输出目录存在
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        
        # 按p值排序并导出
        sorted_df = self.merged_df.sort_values(by='p_val', ascending=True)
        
        # 导出文件
        if out_file.endswith('.xlsx'):
            sorted_df.to_excel(out_file, index=False)
        else:
            sorted_df.to_csv(out_file, sep='\t', index=False)
        
        print(f"注释文件已保存到: {out_file}")

if __name__ == "__main__":
    # 输入输出文件路径
    marker_file = r"D:\nextcloud\pd论文\data\cluster-marker\cluster0_markers.tsv"
    annotation_file = r"D:\nextcloud\pd论文\data\g2mysijfp52-3\Spis_gene_annotation.tsv"  # 基因注释文件
    output_file = r"D:\nextcloud\pd论文\result\annotation1\cluster0_markers_annotated.tsv"
    
    # 处理注释
    anno = AnnotationCellType(marker_file, annotation_file)
    anno.read_marker_file()
    anno.read_annotation_file()
    anno.merge_data()
    anno.export(output_file)























