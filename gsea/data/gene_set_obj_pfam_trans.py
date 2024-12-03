import os
from collections import OrderedDict
import pandas as pd


class Gmt_stat():
    def __init__(self, file_name:str, sepr='\t', ele_sepr=':', _lower=True):
        self.file_name = file_name
        self.sepr = sepr
        self.ele_sepr = ele_sepr
        self.term_number = {}
        self.term_gene_dic = {}
        self.term_name = {}
        self.term_class = {}
        self.head_list = []
        self.unique_genes = set()
        self.get_gmtobj()

    def read_as_pd(self):
        tsv_data = pd.read_csv(self.file_name, sep='\t', header=0, index_col=0)

    def get_gmtobj(self):
        uniq_gene = []
        with open(self.file_name, 'r') as fgmt:
            self.head_list = fgmt.readline().strip('\n').split(self.sepr)
            for line in fgmt.readlines():
                Gene_id, Pfam_number, *Pfam_id_Pfam_description = line.strip('\n').split(self.sepr)
                Gene_id = Gene_id.lower()
                uniq_gene.append(Gene_id)
                for each_pfam in Pfam_id_Pfam_description:
                    pfam_list = self.get_go_info(each_pfam)
                    if not pfam_list:
                        continue
                    if pfam_list[0] not in self.term_name:
                        self.term_name[pfam_list[0]] = pfam_list[1]
                        self.term_gene_dic[pfam_list[0]] = [Gene_id]
                    else:
                        self.term_gene_dic[pfam_list[0]].append(Gene_id)
        self.unique_genes = set(uniq_gene)

    def get_go_info(self, each_pfam):
        try:
            Pfam_id, Pfam_description = each_pfam.split(self.ele_sepr)
        except:
            return False
        return [Pfam_id, Pfam_description]

    def low_case(self, data_list: list):
        return [data.lower() for data in data_list]

    @staticmethod
    def from_gmt(gmt_file: str) -> 'Gmt_stat':
        """从GMT文件创建Gmt_stat对象
        
        GMT文件格式:
        <Pfam_id>\t<Pfam_description>\t<gene1>\t<gene2>...\t<geneN>
        """
        gmt_obj = Gmt_stat.__new__(Gmt_stat)
        gmt_obj.file_name = gmt_file
        gmt_obj.sepr = '\t'
        gmt_obj.ele_sepr = ':'
        gmt_obj.term_number = {}
        gmt_obj.term_gene_dic = {}
        gmt_obj.term_name = {}
        gmt_obj.head_list = []
        gmt_obj.unique_genes = set()

        with open(gmt_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 3:  # 至少需要Pfam_id、Pfam_description和一个基因
                    continue
                    
                pfam_id, pfam_desc, *genes = fields
                genes = [gene.lower() for gene in genes]  # 保持基因名小写
                
                gmt_obj.term_name[pfam_id] = pfam_desc
                gmt_obj.term_gene_dic[pfam_id] = genes
                gmt_obj.unique_genes.update(genes)
                
        return gmt_obj

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

Gene_id	Pfam_number	Pfam_id:Pfam_description

'''

if __name__ == '__main__':
    file_fp = r"E:\coral_experiment\GSEA\data\OA1_4MF.pfam.anno.xls"
    out_fp = r"E:\coral_experiment\GSEA\data\pfam.gmt"
    gmt_obj = Gmt_stat(file_name=file_fp)
    gmt_obj.get_gmtobj()
    exper = ExportGmt(gmt_obj, out_fp)
    exper.export()



