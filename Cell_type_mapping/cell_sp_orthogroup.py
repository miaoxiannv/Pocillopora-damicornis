'''
orthologous_pairs.txt 保存的是同源基因对， 文件格式如下：

Nvec_v1g150206	Spis_XP_022808517_1
Nvec_v1g231956	Spis_XP_022808517_1
Nvec_v1g150206	Fspp_ffun1_m4_18844_m1
Nvec_v1g231956	Fspp_ffun1_m4_18844_m1
'''

import pandas as pd
class Cell_sp_orthogroup:
    '''
    该类用于处理细胞类型与sp的orthologous pairs文件,结果保存为一个字典，字典的键为基因名，字典的值为基因名的列表
    '''
    def __init__(self, orthologous_pairs_file):
        self.orthologous_pairs_file = orthologous_pairs_file
        self.orthologous_pairs_df = pd.read_csv(orthologous_pairs_file, header=None, index_col=False, sep='\t')
        self.orthologous_pairs_df.columns = ['first_gene', 'second_gene']
        self.gene_pair_dict = {}
        self.get_cell_sp_orthogroup()

    def get_cell_sp_orthogroup(self):
        for idx, row in self.orthologous_pairs_df.iterrows():
            #如果字典中没有该基因对，则添加该基因对，如果字典中有该基因对，则将第二个基因添加到列表中
            if row['first_gene'] not in self.gene_pair_dict:
                self.gene_pair_dict[row['first_gene']] = [row['second_gene']]
            else:
                self.gene_pair_dict[row['first_gene']].append(row['second_gene'])
            #以第二个基因作为键，第一个基因作为值，添加到字典中
            if row['second_gene'] not in self.gene_pair_dict:
                self.gene_pair_dict[row['second_gene']] = [row['first_gene']]
            else:
                self.gene_pair_dict[row['second_gene']].append(row['first_gene'])




