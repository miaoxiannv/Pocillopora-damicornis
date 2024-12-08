import os
from collections import OrderedDict, defaultdict
import pandas as pd


class Gmt_stat():
    def __init__(self, file_name: str, sepr='\t', ele_sepr=';', _lower=True):
        self.file_name = file_name
        self.sepr = sepr
        self.ele_sepr = ele_sepr
        self.term_number = {}
        self.term_gene_dic = {}
        self.term_name = {}
        self.head_list = []
        self.unique_genes = set()
        self.get_gmtobj()

    def read_as_pd(self):
        tsv_data = pd.read_csv(self.file_name, sep='\t', header=0, index_col=0)

    def get_gmtobj(self):
        with open(self.file_name, 'r') as fgmt:
            self.head_list = fgmt.readline().strip('\n').split(self.sepr)
            for line in fgmt.readlines():
                term_id, term_name, term_number, *gene_list = line.strip('\n').split(self.sepr)
                term_id = term_id.upper()
                self.term_name[term_id] = term_name
                self.term_number[term_id] = term_number

                gene_list = self.low_case(list(gene_list))
                self.term_gene_dic[term_id] = list(gene_list)

                self.unique_genes.update(set(gene_list))

    def low_case(self, data_list: list):
        return [data.lower() for data in data_list]

    @staticmethod
    def from_gmt(gmt_file: str) -> 'Gmt_stat':
        """从GMT文件创建Gmt_stat对象
        
        GMT文件格式:
        <term>\t<term_name>\t<gene1>\t<gene2>...\t<geneN>
        """
        gmt_obj = Gmt_stat.__new__(Gmt_stat)
        gmt_obj.file_name = gmt_file
        gmt_obj.sepr = '\t'
        gmt_obj.ele_sepr = ';'
        gmt_obj.term_number = {}
        gmt_obj.term_gene_dic = {}
        gmt_obj.term_name = {}
        gmt_obj.head_list = []
        gmt_obj.unique_genes = set()

        with open(gmt_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 3:  # 至少需要term、term_name和一个基因
                    continue
                    
                term_id, term_name, *genes = fields
                term_id = term_id.upper()
                
                gmt_obj.term_name[term_id] = term_name
                gmt_obj.term_gene_dic[term_id] = genes
                gmt_obj.term_number[term_id] = str(len(genes))
                gmt_obj.unique_genes.update(genes)
                
        return gmt_obj


class Gmt_stat_gene():
    def __init__(self, file_name: str, sepr='\t', ele_sepr=';', _lower=True, term_sepr=" | "):
        self.file_name = file_name
        self.sepr = sepr
        self.ele_sepr = ele_sepr
        self.term_sepr = term_sepr
        self.term_number = {}
        self.term_gene_dic = defaultdict(set)
        self.term_name = {}
        self.head_list = []
        self.unique_genes = set()
        self.get_gmtobj()

    def get_gmtobj(self):
        with open(self.file_name, 'r') as fgmt:
            self.head_list = fgmt.readline().strip('\n').split(self.sepr)
            for line in fgmt.readlines():
                Query_id, Subject_id, KO_ID, KO_NAME, KO_DEFINITION, KO_EC, KO_PATHWAY = line.strip('\n').split(self.sepr)
                gene = Query_id.lower()
                self.unique_genes.update([gene])
                if KO_PATHWAY == "-":
                    continue
                each_KO_PATHWAY = KO_PATHWAY.split(self.term_sepr)
                for kp in each_KO_PATHWAY:
                    res = self.get_kegg_info(kp)
                    if isinstance(res, list):
                        kegg_id = res.pop(0).lower()
                        kegg_name = self.ele_sepr.join(res)
                        self.term_gene_dic[kegg_id].update([gene])
                        self.term_name[kegg_id] = kegg_name
        self.term_gene_dic = {k: list(v) for k, v in self.term_gene_dic.items()}

    def get_kegg_info(self, each_pfam):
        try:
            kegg_id, kegg_class, kegg_subclass, kegg_term= each_pfam.split(self.ele_sepr)
        except:
            return False
        return [kegg_id, kegg_class, kegg_subclass, kegg_term]

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
KO_Pathway_Level3	Pathway_Name	Reads_Num	Reads_IDs.


'''
if __name__ == '__main__':
    file_fp = r"D:\nextcloud\pd论文\data\test\gsea\P2.KEGG.filter.m8.anno.xls"
    out_fp = r"D:\nextcloud\pd论文\data\test\gsea\kegg.gmt"
    gmt_obj = Gmt_stat(file_name=file_fp)
    gmt_obj.get_gmtobj()
    exper = ExportGmt(gmt_obj, out_fp)
    exper.export()
