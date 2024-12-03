import os, sys

'''
输入的文件分为两列，第一列是细胞类型，第二列是marker基因
文件名可能是csv或者是tsv格式
'''

class ReadMarker:
    def __init__(self, filename):
        self.filename = filename
        self.marker_dict = {}
        self.read_file()

    def read_file(self):
        if self.filename.endswith('.csv'):
            self.read_csv()
        elif self.filename.endswith('.tsv'):
            self.read_tsv()
        else:
            print('文件格式错误，只支持csv或者tsv格式文件')
            sys.exit(1)

    def read_csv(self):
        with open(self.filename, 'r') as f:
            for line in f:
                line = line.strip().lower()
                # replace "/" and "-" with "_"
                line = line.replace("/", "_").replace("-", "_")
                if not line:
                    continue
                line = line.split(',')
                if len(line) != 2:
                    print('文件格式错误，csv文件应该有两列')
                    sys.exit(1)
                cell_type = line[0].strip()
                marker_gene = line[1].strip()
                if cell_type not in self.marker_dict:
                    self.marker_dict[cell_type] = []
                self.marker_dict[cell_type].append(marker_gene)

            self.marker_dict = {cell_type: set(genes) for cell_type, genes in self.marker_dict.items()}

    def read_tsv(self):
        with open(self.filename, 'r') as f:
            for line in f:
                line = line.strip().lower()
                # replace "/" and "-" with "_"
                line = line.replace("/", "_").replace("-", "_")
                if not line:
                    continue
                line = line.split('\t')
                if len(line) != 2:
                    print('文件格式错误，tsv文件应该有两列')
                    sys.exit(1)
                cell_type = line[0].strip()
                marker_gene = line[1].strip()
                if cell_type not in self.marker_dict:
                    self.marker_dict[cell_type] = set()
                self.marker_dict[cell_type].update({marker_gene})

            #self.marker_dict = {cell_type: set(genes) for cell_type, genes in self.marker_dict.items()}

    def get_marker(self):
        return self.marker_dict
    

    

