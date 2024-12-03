import os
from collections import OrderedDict
import pandas as pd


class Gmt_stat():
    def __init__(self, file_name:str, sepr='\t', _lower=True):
        self.file_name = file_name
        self.sepr = sepr
        self.term_number = {}
        self.term_gene_dic = {}
        self.term_name = {}
        self.term_class = {}
        self.head_list = []
        self.unique_genes = set()
        self.get_gmtobj()
        self.ele_sepr = ";"

    def read_as_pd(self):
        tsv_data = pd.read_csv(self.file_name, sep='\t', header=0, index_col=0)

    def get_gmtobj(self):
        uniq_gene = []
        with open(self.file_name, 'r') as fgmt:
            self.head_list = fgmt.readline().strip('\n').split(self.sepr)
            for line in fgmt.readlines():
                Gene_id, Identity, E_value, KOG_gene_id, KOG_num, Functional_description, Functional_class, Class_description = line.strip('\n').split(self.sepr)
                Gene_id = Gene_id.lower()
                uniq_gene.append(Gene_id)
                if KOG_num not in self.term_name:
                    self.term_name[KOG_num] = Functional_description
                    self.term_gene_dic[KOG_num] = [Gene_id]
                    self.term_class[KOG_num] = Class_description
                else:
                    self.term_gene_dic[KOG_num].append(Gene_id)
        self.unique_genes = set(uniq_gene)

    @staticmethod
    def from_gmt(gmt_file: str) -> 'Gmt_stat':
        """从GMT文件创建Gmt_stat对象
        
        GMT文件格式:
        <KOG_num>\t<Functional_description>\t<gene1>\t<gene2>...\t<geneN>
        """
        gmt_obj = Gmt_stat.__new__(Gmt_stat)
        gmt_obj.file_name = gmt_file
        gmt_obj.sepr = '\t'
        gmt_obj.ele_sepr = ';'
        gmt_obj.term_number = {}
        gmt_obj.term_gene_dic = {}
        gmt_obj.term_name = {}
        gmt_obj.term_class = {}
        gmt_obj.head_list = []
        gmt_obj.unique_genes = set()

        with open(gmt_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 3:  # 至少需要KOG_num、Functional_description和一个基因
                    continue
                    
                kog_num, func_desc, *genes = fields
                genes = [gene.lower() for gene in genes]  # 保持基因名小写
                
                gmt_obj.term_name[kog_num] = func_desc
                gmt_obj.term_gene_dic[kog_num] = genes
                gmt_obj.term_class[kog_num] = "KOG"  # KOG分类
                gmt_obj.unique_genes.update(genes)
                
        return gmt_obj

    def get_go_info(self, go_3):
        try:
            go_id, go_description, go_class = go_3.split(self.ele_sepr)
        except:
            return False
        return [go_id, go_description, go_class]

    def low_case(self, data_list: list):
        return [data.lower() for data in data_list]

class ExportGmt(object):
    def __init__(self, gmt_obj:Gmt_stat, out_fp:str):
        self.gmt_obj = gmt_obj
        self.out_fp = out_fp

    def export(self):
        with open(self.out_fp, 'w') as fw:
            for term, genes in self.gmt_obj.term_gene_dic.items():
                term_ = "\t".join([term, gmt_obj.term_name[term]])
                genes = "\t".join(genes)
                fw.write("\t".join([term_, genes]))
                fw.write("\n")



'''

Gene_id	Identity	E_value	KOG_gene_id	KOG_num	Functional_description	Functional_class	Class_description

'''
if __name__ == '__main__':
    file_fp = r"E:\coral_experiment\GSEA\data\OA1_4MF.KOG.filter.m8.anno.xls"
    out_fp = r"E:\coral_experiment\GSEA\data\kog.gmt"
    gmt_obj = Gmt_stat(file_name=file_fp)
    gmt_obj.get_gmtobj()
    exper = ExportGmt(gmt_obj, out_fp)
    exper.export()



