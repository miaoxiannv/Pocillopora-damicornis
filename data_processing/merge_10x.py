'''
this script is used to merge 10x data from different species
input directory: barcode.tsv, feature.tsv, matrix.mtx
output directory: barcodes.tsv, features.tsv, matrix.mtx
'''

import gzip
import os

'''
barcode.tsv 文件内容如下：
AAACCTGCAGCCTGTG-1
AAACCTGGTCCGTCAG-1
AAACCTGGTTCTGTTT-1
AAACCTGTCCGCTGTT-1
AAACCTGTCTGGTATG-1
AAACGGGAGTCGAGTG-1

'''


class Barcode:
    def __init__(self, file_path):
        self.file_path = file_path
        self.barcode = self.read_barcode()

    def read_barcode(self):
        if self.file_path.endswith('.gz'):
            with gzip.open(self.file_path, 'rt', encoding='utf-8') as f:
                barcode = [line.strip() for line in f]
        else:
            with open(self.file_path, 'r', encoding='utf-8') as f:
                barcode = [line.strip() for line in f]
        return barcode


'''
feature.tsv 文件内容如下：
transcript_HQ_OA1_1SP_transcript10903/f12p0/2901	transcript_HQ_OA1_1SP_transcript10903/f12p0/2901	Gene Expression
transcript_HQ_OA1_1SP_transcript12719/f24p0/2754	transcript_HQ_OA1_1SP_transcript12719/f24p0/2754	Gene Expression
transcript_HQ_OA1_1SP_transcript14531/f6p0/2646	transcript_HQ_OA1_1SP_transcript14531/f6p0/2646	Gene Expression
transcript_HQ_OA1_1SP_transcript13629/f3p0/2697	

'''


class Feature:
    def __init__(self, file_path):
        self.file_path = file_path
        self.feature = self.read_feature()

    def read_feature(self):
        if self.file_path.endswith('.gz'):
            with gzip.open(self.file_path, 'rt', encoding='utf-8') as f:
                feature = [line.strip() for line in f]
        else:
            with open(self.file_path, 'r', encoding='utf-8') as f:
                feature = [line.strip() for line in f]
        return feature


'''
%%MatrixMarket matrix coordinate integer general
%metadata_json: {"software_version": "cellranger-8.0.1", "format_version": 2}
22408 2122 156897
53 1 2
1103 1 1
1426 1 1
1652 1 1

其中，第一行和第二行是注释行，第三行是矩阵的维度信息，第四行开始是矩阵的数据，每一行的数据格式是：行号 列号 数据值
数据有三列，分别是feature_id, barcode_id, count
'''


class Matrix:
    def __init__(self, file_path, sepr=' '):
        self.file_path = file_path
        self.sepr = sepr
        self.matrix = self.read_matrix()

    def read_matrix(self):
        if self.file_path.endswith('.gz'):
            with gzip.open(self.file_path, 'rt', encoding='utf-8') as f:
                lines = f.readlines()[3:]  # 跳过前三行
                matrix = [line.strip().split(self.sepr) for line in lines]
        else:
            with open(self.file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()[3:]  # 跳过前三行
                matrix = [line.strip().split(self.sepr) for line in lines]
        return matrix


# 将数据保存为一个字典，结构如下：{barcode_id: {feature_id: count}}
class MergeData:
    def __init__(self, barcode, feature, matrix):
        self.barcode = barcode
        self.feature = feature
        self.matrix = matrix
        self.data = self.trans_data()

    def trans_data(self):
        data = {}
        for line in self.matrix.matrix:
            barcode_id = line[1]
            feature_id = line[0]
            # 获取barcode_id对应的barcode name
            barcode_name = self.barcode.barcode[int(barcode_id) - 1]
            #   获取feature_id对应的feature name
            feature_name = self.feature.feature[int(feature_id) - 1]
            count = line[2]
            # 如果barcode_id不在data中，则初始化一个空字典
            if barcode_name not in data:
                data[barcode_name] = {}
            # 如果feature_id不在data[barcode_id]中，则初始化为0
            if feature_name not in data[barcode_name]:
                data[barcode_name][feature_name] = 0
            data[barcode_name][feature_name] += int(count)
        return data

    # 将一个新的MergeData对象与当前的MergeData对象合并
    def merge_data(self, new_data):
        for barcode in new_data.data:
            if barcode not in self.data:
                self.data[barcode] = new_data.data[barcode]
            else:
                for feature in new_data.data[barcode]:
                    if feature not in self.data[barcode]:
                        self.data[barcode][feature] = new_data.data[barcode][feature]
                    else:
                        self.data[barcode][feature] += new_data.data[barcode][feature]
        return self.data


# 将MergeData.data保存为为三个文件，分别是：matrix.mtx, features.tsv, barcodes.tsv
class SaveData:
    def __init__(self, data, file_path='D:\Bioinformatics\数据整合\最终版数据\下游分析\\test001'):
        self.data = data
        self.file_path = file_path
        # barcode file path, feature file path, matrix file path
        self.barcode_file = os.path.join(self.file_path, 'barcodes.tsv')
        self.feature_file = os.path.join(self.file_path, 'features.tsv')
        self.matrix_file = os.path.join(self.file_path, 'matrix.mtx')
        self.save_data()

    # self.data为MergeData.data，遍历self.data，将数据保存到三个文件中
    def save_data(self):
        # 获得所有的barcode name
        barcode = list(self.data.keys())
        barcode_to_id = {barcode: i + 1 for i, barcode in enumerate(barcode)}
        # 获得所有的feature name
        feature = list(set([feature for barcode in self.data for feature in self.data[barcode]]))
        feature_to_id = {feature: i + 1 for i, feature in enumerate(feature)}

        # 将barcode name保存到文件中
        with open(self.barcode_file, 'w') as f:
            for bc in barcode:
                f.write(bc + '\n')
        # 将feature name保存到文件中
        with open(self.feature_file, 'w') as f:
            for ft in feature:
                f.write(ft + '\n')
        # 将数据保存到matrix.mtx文件中, 格式为：feature_id barcode_id count
        with open(self.matrix_file, 'w') as f:
            content = []
            for bc in barcode:
                for ft in feature:
                    if ft in self.data[bc]:
                        # feature name对应的feature id, barcode name对应的barcode id, count
                        content.append([feature_to_id[ft], barcode_to_id[bc], self.data[bc][ft]])

            f.write('%%MatrixMarket matrix coordinate integer general\n')
            f.write('%metadata_json: {"software_version": "cellranger-8.0.1", "format_version": 2}\n')
            f.write(f'{len(feature)} {len(barcode)} {len(content)}\n')
            for line in content:
                f.write(' '.join([str(i) for i in line]) + '\n')


if __name__ == '__main__':
    # 定义5个物种数据集的基础路径
    data_paths = {
        'species1': {
            'barcode': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\SP-data\\barcodes.tsv',
            'feature': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\SP-data\\features.tsv',
            'matrix': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\SP-data\matrix.mtx'
        },
        'species2': {
            'barcode': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\PD-data\\barcodes.tsv',
            'feature': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\PD-data\\features.tsv',
            'matrix': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\PD-data\matrix.mtx'
        },
        'species3': {
            'barcode': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\PV-data\\barcodes.tsv',
            'feature': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\PV-data\\features.tsv',
            'matrix': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\PV-data\\matrix.mtx'
        },
        'species4': {
            'barcode': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\FL-data\\barcodes.tsv',
            'feature': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\FL-data\\features.tsv',
            'matrix': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\FL-data\\matrix.mtx'
        },
        'species5': {
            'barcode': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\SD-data\\barcodes.tsv',
            'feature': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\SD-data\\features.tsv',
            'matrix': 'D:\Bioinformatics\数据整合\最终版数据\下游分析\SD-data\\matrix.mtx'
        }

    }

    # 初始化数据列表
    merge_data_list = []

    # 读取所有数据集
    for species, paths in data_paths.items():
        try:
            print(f"Processing {species}...")
            barcode = Barcode(paths['barcode'])
            feature = Feature(paths['feature'])
            matrix = Matrix(paths['matrix'])

            # 创建MergeData对象
            data = MergeData(barcode, feature, matrix)
            merge_data_list.append(data)
            print(f"Successfully processed {species}")
        except Exception as e:
            print(f"Error processing {species}: {str(e)}")
            continue

    if not merge_data_list:
        print("No data was successfully processed. Exiting...")
        exit()

    # 合并所有数据集
    final_data = merge_data_list[0]
    for data in merge_data_list[1:]:
        final_data.merge_data(data)

    # 保存合并后的数据
    out_file_path = r'D:\Bioinformatics\数据整合\最终版数据\下游分析\test001'
    os.makedirs(out_file_path, exist_ok=True)  # 确保输出目录存在
    save_data = SaveData(final_data.data, out_file_path)

    print('Merging completed successfully')
