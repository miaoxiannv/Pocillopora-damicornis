from typing import Union

'''
p_val	avg_log2FC	pct.1	pct.2	p_val_adj	cluster	gene
0.000359921	12.02402076	0.011	0	1	0	LOC113670121
0.000516433	12.02238502	0.01	0	1	0	LOC131799609
0.000741732	11.85736729	0.01	0	1	0	transcript-HQ-P2-transcript20480/f2p0/587
'''

import pandas as pd
import numpy as np

class MarkerGeneReader:
    """读取marker gene文件的类
    
    Attributes:
        file_path (str): marker gene文件路径
        df (pd.DataFrame): 存储marker gene数据的DataFrame
        clusters (list): 所有cluster的列表
        genes (set): 所有基因的集合
    """
    
    def __init__(self, file_path: Union[str, pd.DataFrame], 
                 p_val_threshold: float = 0.05,
                 log2fc_threshold: float = 1.0):
        """初始化MarkerGeneReader
        
        Args:
            file_path (str): marker gene文件路径
            p_val_threshold (float): p值阈值，默认0.05
            log2fc_threshold (float): log2FC阈值，默认1.0
        """
        self.file_path = file_path
        self.p_val_threshold = p_val_threshold
        self.log2fc_threshold = log2fc_threshold
        
        # 读取数据
        if isinstance(file_path, str):
            self.df = pd.read_csv(file_path, sep='\t')
        elif isinstance(file_path, pd.DataFrame):
            self.df = file_path
        else:
            raise ValueError("file_path must be a string or a pandas DataFrame")
        
        # 获取基本信息
        self.clusters = sorted(self.df['cluster'].unique())
        self.genes = set(self.df['gene'])
        
    def get_cluster_markers(self, cluster_id: int, filtered: bool = True) -> pd.DataFrame:
        """获取特定cluster的marker基因
        
        Args:
            cluster_id (int): cluster的ID
            filtered (bool): 是否根据p值和log2FC阈值过滤，默认True
            
        Returns:
            pd.DataFrame: 该cluster的marker基因数据
        """
        cluster_data = self.df[self.df['cluster'] == cluster_id]
        
        if filtered:
            cluster_data = cluster_data[
                (cluster_data['p_val'] < self.p_val_threshold) & 
                (cluster_data['avg_log2FC'].abs() > self.log2fc_threshold)
            ]
        
        return cluster_data.sort_values('p_val')
    # get negative cluster markers
    def get_negative_cluster_markers(self, cluster_id: int, filtered: bool = True) -> pd.DataFrame:
        cluster_data = self.df[self.df['cluster'] == cluster_id]
        if filtered:
            cluster_data = cluster_data[
                (cluster_data['p_val'] < self.p_val_threshold) & 
                (cluster_data['avg_log2FC'] < -self.log2fc_threshold)
            ]
        return cluster_data.sort_values('p_val', ascending=True)

    
    def get_gene_info(self, gene_id: str) -> pd.Series:
        """获取特定基因的信息
        
        Args:
            gene_id (str): 基因ID
            
        Returns:
            pd.Series: 该基因的所有信息，如果不存在返回None
        """
        gene_data = self.df[self.df['gene'] == gene_id]
        return gene_data.iloc[0] if not gene_data.empty else None
    
    def get_significant_markers(self) -> pd.DataFrame:
        """获取所有显著的marker基因
        
        Returns:
            pd.DataFrame: 满足p值和log2FC阈值的marker基因
        """
        return self.df[
            (self.df['p_val'] < self.p_val_threshold) & 
            (self.df['avg_log2FC'].abs() > self.log2fc_threshold)
        ].sort_values(['cluster', 'p_val'])
    
    def get_stats(self) -> dict:
        """获取基本统计信息
        
        Returns:
            dict: 包含基本统计信息的字典
        """
        sig_markers = self.get_significant_markers()
        
        stats = {
            'total_genes': len(self.genes),
            'total_clusters': len(self.clusters),
            'significant_markers': len(sig_markers),
            'markers_per_cluster': sig_markers.groupby('cluster').size().to_dict()
        }
        
        return stats
    
    def to_dict(self) -> dict:
        """将marker基因数据转换为字典格式
        
        Returns:
            dict: {cluster_id: {gene_id: marker_info}}的嵌套字典
        """
        marker_dict = {}
        for cluster in self.clusters:
            cluster_data = self.get_cluster_markers(cluster)
            marker_dict[cluster] = cluster_data.set_index('gene').to_dict('index')
        
        return marker_dict
    
    def get_gene_list_from_df(self, df: pd.DataFrame) -> list:
        """从DataFrame中获取基因列表
        
        Args:
            df (pd.DataFrame): 包含基因信息的DataFrame
            
        Returns:
            list: 基因列表
        """
        return df['gene'].tolist()
    
    def save_gene_list(self, gene_list: list, file_path: str):
        """保存基因列表到文件
        
        Args:
            gene_list (list): 基因列表
            file_path (str): 保存路径
        """
        with open(file_path, 'w') as f:
            for gene in gene_list:
                f.write(f"{gene}\n")


if __name__ == '__main__':
    # 创建测试数据
    test_data = {
        'p_val': [0.001, 0.01, 0.005, 0.02],
        'gene': ['gene1', 'gene2', 'gene3', 'gene4'],
        'cluster': [0, 0, 1, 1],
        'avg_log2FC': [2.5, 1.8, 3.2, 2.1],
        'p_val_adj': [0.001, 0.01, 0.005, 0.02],
        'pct.1': [0.8, 0.7, 0.9, 0.6],
        'pct.2': [0.2, 0.3, 0.1, 0.2]
    }
    test_df = pd.DataFrame(test_data)
    
    # 将 MarkerGenes 改为 MarkerGeneReader
    marker_genes = MarkerGeneReader(test_df)
    print(marker_genes.get_cluster_markers(0))
    print(marker_genes.get_negative_cluster_markers(0))
    print(marker_genes.get_gene_info('gene1'))
    print(marker_genes.get_significant_markers())
    print(marker_genes.get_stats())
    print(marker_genes.to_dict())
    print(marker_genes.get_gene_list_from_df(test_df))
    marker_genes.save_gene_list(marker_genes.get_gene_list_from_df(test_df), 'test_genes.txt')







