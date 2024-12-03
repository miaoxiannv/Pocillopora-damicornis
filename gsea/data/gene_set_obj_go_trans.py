
'''
go annotation xlsx file format:
Gene_id	GO_number	GO_id; GO_description; GO_class	
transcript_HQ_P2_transcript0/f8p0/9897	1	GO:0005515; protein binding; molecular_function	
transcript_HQ_P2_transcript1/f2p0/9329	10	GO:0006914; autophagy; biological_process	GO:0005198; NA
transcript_HQ_P2_transcript10/f2p0/8399	4	GO:0004222; metalloendopeptidase activity; molecular_function	GO:0003677; DNA binding; molecular_function
transcript_HQ_P2_transcript100/f5p0/6808	8	GO:0016192; vesicle-mediated transport; biological_process	GO:0006886; intracellular protein transport; biological_process
transcript_HQ_P2_transcript1000/f3p0/4849	7	GO:0016787; hydrolase activity; molecular_function	GO:0003677; DNA binding; molecular_function
transcript_HQ_P2_transcript10000/f12p0/2537	4	GO:0005666; DNA-directed RNA polymerase III complex; cellular_component	GO:0006383; transcription from RNA polymerase III promoter; biological_process
transcript_HQ_P2_transcript10001/f2p0/2514	4	GO:0015031; protein transport; biological_process	GO:0035091; phosphatidylinositol binding; cellular_component

'''



import os
from collections import OrderedDict
import pandas as pd


class Gmt_stat():
    def __init__(self, file_name:str, sepr='\t', ele_sepr=';', _lower=True):
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
        # rename the term_class keys
        #self.term_abbr = {"molecular_function": "mf",
         #                 "biological_process": "bp",
          #                "cellular_component": "cc"}


    def read_as_pd(self):
        tsv_data = pd.read_csv(self.file_name, sep='\t', header=0, index_col=0)

    def get_gmtobj(self):
        uniq_gene = []
        with open(self.file_name, 'r') as fgmt:
            self.head_list = fgmt.readline().strip('\n').split(self.sepr)
            for line in fgmt.readlines():
                Gene_id, GO_number, *GO_id_GO_description_GO_class = line.strip('\n').split(self.sepr)
                Gene_id = Gene_id.lower()
                uniq_gene.append(Gene_id)
                for each_go in GO_id_GO_description_GO_class:
                    go_list = self.get_go_info(each_go)
                    if not go_list:
                        continue
                    if go_list[0] not in self.term_name:
                        self.term_name[go_list[0]] = go_list[1]
                        self.term_class[go_list[0]] = go_list[2]
                        self.term_gene_dic[go_list[0]] = [Gene_id]
                    else:
                        self.term_gene_dic[go_list[0]].append(Gene_id)
        self.unique_genes = set(uniq_gene)

    def get_go_info(self, go_3):
        try:
            go_id, go_description, go_class = go_3.strip().split(self.ele_sepr)
        except:
            return False
        return [go_id.strip(), go_description.strip(), go_class.strip()]

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
        gmt_obj.term_class = {}
        gmt_obj.head_list = []
        gmt_obj.unique_genes = set()

        with open(gmt_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 3:  # 至少需要term、term_name和一个基因
                    continue
                    
                term, term_name = fields[0:2]
                genes = fields[2:]
                
                gmt_obj.term_name[term] = term_name
                gmt_obj.term_gene_dic[term] = genes
                gmt_obj.unique_genes.update(genes)
                # 由于GMT文件通常不包含GO分类信息，默认设为biological_process
                gmt_obj.term_class[term] = "biological_process"
                
        return gmt_obj

    def split_by_term_class(self) -> dict:
        """将当前Gmt_stat对象按term_class分割成三个子对象
        
        Returns:
            dict: 包含三个子Gmt_stat对象的字典
                {
                    'bp': Gmt_stat object (biological_process),
                    'mf': Gmt_stat object (molecular_function),
                    'cc': Gmt_stat object (cellular_component)
                }
        """
        # 初始化三个子对象
        sub_objs = {}
        for class_type in ['bp', 'mf', 'cc']:
            sub_obj = Gmt_stat.__new__(Gmt_stat)
            sub_obj.file_name = self.file_name
            sub_obj.sepr = self.sepr
            sub_obj.ele_sepr = self.ele_sepr
            sub_obj.term_number = {}
            sub_obj.term_gene_dic = {}
            sub_obj.term_name = {}
            sub_obj.term_class = {}
            sub_obj.head_list = self.head_list
            sub_obj.unique_genes = set()
            sub_objs[class_type] = sub_obj
        
        # 映射完整类名到简写
        class_map = {
            "biological_process": "bp",
            "molecular_function": "mf",
            "cellular_component": "cc"
        }
        
        # 根据term_class分配数据
        for term, class_type in self.term_class.items():
            if class_type not in class_map:
                continue
                
            class_abbr = class_map[class_type]
            sub_obj = sub_objs[class_abbr]
            
            # 复制相关数据到子对象
            if term in self.term_gene_dic:
                sub_obj.term_gene_dic[term] = self.term_gene_dic[term]
                sub_obj.unique_genes.update(self.term_gene_dic[term])
            
            if term in self.term_name:
                sub_obj.term_name[term] = self.term_name[term]
                
            if term in self.term_number:
                sub_obj.term_number[term] = self.term_number[term]
                
            sub_obj.term_class[term] = class_type
            
        return sub_objs




class ExportGmt(object):
    def __init__(self, gmt_obj:Gmt_stat, out_fp:str):
        self.gmt_obj = gmt_obj
        self.out_fp = out_fp

    def export(self):
        # 创建输出目录
        out_dir = os.path.dirname(self.out_fp)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            
        with open(self.out_fp + "_bp.gmt", 'w') as fbp:
            with open(self.out_fp + "_cc.gmt", 'w') as fcc:
                with open(self.out_fp + "_mf.gmt", 'w') as fmf:
                    fp_map = {"molecular_function": fmf,
                              "cellular_component": fcc,
                              "biological_process": fbp}
                    for term, genes in self.gmt_obj.term_gene_dic.items():
                        term_ = "\t".join([term, gmt_obj.term_name[term]])
                        genes = "\t".join(genes)
                        fw = fp_map[gmt_obj.term_class[term]]
                        fw.write("\t".join([term_, genes]))
                        fw.write("\n")


'''
Gene_id	GO_number	GO_id; GO_description; GO_class

'''
if __name__ == '__main__':
    file_fp = r"D:\Bioinformatics\数据整合\GSEA\GSEA\new-data\GO.anno.xls"
    out_fp = r"D:\Bioinformatics\数据整合\GSEA\result\GO\go.gmt"
    gmt_obj = Gmt_stat(file_name=file_fp)
    gmt_obj.get_gmtobj()
    sub_objs = gmt_obj.split_by_term_class()
    for class_type, sub_obj in sub_objs.items():
        out_fp_class = out_fp + "_" + class_type + ".gmt"
        exper = ExportGmt(sub_obj, out_fp_class)
        exper.export()
    exper = ExportGmt(gmt_obj, out_fp)
    exper.export()


