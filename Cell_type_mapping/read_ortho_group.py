import sys

'''
输入文件分三列，第一列是Orthogroup ID，第二列第一个物种对应的基因ID的集合，以","分割，第三列第二个物种对应的基因ID的集合，以","分割
列与列之间以"\t"分割
文件含有表头：
Orthogroup	FL	SP-cell
OG0000000	transcript_HQ_P2_transcript1877_f2p0_4197_1, transcript_HQ_P2_transcript2511_f2p0_3915_1, transcript_HQ_P2_transcript4149_f2p0_3397_1	Spis10168_1, Spis10297_1, Spis10529_1, Spis10717_1, 
'''

class ReadOrthoGroup:
    def __init__(self, filename, sp1_col, sp2_col):
        try:
            self.filename = filename
            self.ortho_group_dict = {}
            self.gene_mapping_dict_sp1 = {}
            self.gene_mapping_dict_sp2 = {}
            self.sp1_all_genes = set()
            self.sp2_all_genes = set()
            self.sp1_col = sp1_col
            self.sp2_col = sp2_col
            self.read_tsv()
        except Exception as e:
            print(f"初始化过程中出错: {str(e)}")
            sys.exit(1)
    
    def read_tsv(self):
        with open(self.filename, 'r') as f:
            #跳过第一行
            f.readline()

            for line_num, line in enumerate(f, 2):
                line = line.strip('\n')
                if not line:
                    continue
                line = line.split('\t')
                if len(line) != 3:
                    print(f'文件格式错误，tsv文件应该有三列，但在第{line_num}行发现{len(line)}列')
                    print(f'问题行内容: {line}')
                    sys.exit(1)
                ortho_group = line[0].strip()
                gene_id_1 = line[self.sp1_col].strip().lower().split(',')
                gene_id_1 = [gene.strip() for gene in gene_id_1]
                gene_id_2 = line[self.sp2_col].strip().lower().split(',')
                gene_id_2 = [gene.strip() for gene in gene_id_2]

                self.ortho_group_dict[ortho_group]  = [gene_id_1, gene_id_2]
                self.sp1_all_genes.update(gene_id_1)
                self.sp2_all_genes.update(gene_id_2)
                if gene_id_1:
                    for gene in gene_id_1:
                        gene = gene.lower()
                        if gene_id_2 and gene:
                            for gene2 in gene_id_2:
                                if gene2:
                                    gene2 = gene2.lower()
                                    if gene not in self.gene_mapping_dict_sp1:
                                        self.gene_mapping_dict_sp1[gene] = set()
                                    self.gene_mapping_dict_sp1[gene].update({gene2})
                if gene_id_2:
                    for gene in gene_id_2:
                        gene = gene.lower()
                        if gene_id_1 and gene:
                            for gene1 in gene_id_1:
                                if gene1:
                                    gene1 = gene1.lower()
                                    if gene not in self.gene_mapping_dict_sp2:
                                        self.gene_mapping_dict_sp2[gene] = set()
                                    self.gene_mapping_dict_sp2[gene].update({gene1})

    def get_gene_mapping_dict_sp1(self):
        return self.gene_mapping_dict_sp1
    
    def get_gene_mapping_dict_sp2(self):
        return self.gene_mapping_dict_sp2

