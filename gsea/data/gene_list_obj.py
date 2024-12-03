class GeneList_Obj(object):
    def __init__(self, fn, header=False):
        self.header = header
        self.header_info = []
        self.fn = fn
        self.gene_list = []
        self.readfile()
        self.low_case()

    def readfile(self):
        with open(self.fn, 'r') as gf:
            if self.header:
                self.header_info = gf.readline()
            for line in gf.readlines():
                line = line.strip('\n')
                if line:
                    self.gene_list.append(line)

    def low_case(self):
        self.gene_list = [data.lower() for data in self.gene_list]
