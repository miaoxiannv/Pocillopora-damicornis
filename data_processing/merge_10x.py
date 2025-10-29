"""
10X Genomics data merger.
This script merges 10X Genomics data from different species into a unified format.
It processes barcode, feature, and matrix files from multiple species and combines them.

Input files required for each species:
    - barcode.tsv: Contains cell barcodes
    - feature.tsv: Contains gene/feature information
    - matrix.mtx: Contains expression matrix in sparse format

Output files:
    - barcodes.tsv: Combined cell barcodes
    - features.tsv: Combined gene/feature list
    - matrix.mtx: Combined expression matrix

Usage:
    python merge_10x.py

Author: Shengyao Zhang
Date: 2024-12-19
Version: 1.0
"""
import gzip
import os


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


class Matrix:
    def __init__(self, file_path, sepr=' '):
        self.file_path = file_path
        self.sepr = sepr
        self.matrix = self.read_matrix()

    def read_matrix(self):
        if self.file_path.endswith('.gz'):
            with gzip.open(self.file_path, 'rt', encoding='utf-8') as f:
                lines = f.readlines()[3:]
                matrix = [line.strip().split(self.sepr) for line in lines]
        else:
            with open(self.file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()[3:]
                matrix = [line.strip().split(self.sepr) for line in lines]
        return matrix


class MergeData:
    def __init__(self, barcode, feature, matrix):
        self.barcode = barcode
        self.feature = feature
        self.matrix = matrix
        self.data = self.trans_data()

    def trans_data(self):
        data = {}
        for line in self.matrix.matrix:
            try:
                barcode_id = int(line[1])
                feature_id = int(line[0])
                
                barcode_name = self.barcode.barcode[barcode_id - 1]
                feature_name = self.feature.feature[feature_id - 1]
                count = int(line[2])
                
                if barcode_name not in data:
                    data[barcode_name] = {}
                if feature_name not in data[barcode_name]:
                    data[barcode_name][feature_name] = 0
                data[barcode_name][feature_name] += count
                
            except (IndexError, ValueError) as e:
                print(f"Error processing line: {line}")
                print(f"Error details: {str(e)}")
                continue
                
        return data
    
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


class SaveData:
    def __init__(self, data, file_path='./output'):
        self.data = data
        self.file_path = file_path
        self.barcode_file = os.path.join(self.file_path, 'barcodes.tsv')
        self.feature_file = os.path.join(self.file_path, 'features.tsv')
        self.matrix_file = os.path.join(self.file_path, 'matrix.mtx')
        self.save_data()
    
    def save_data(self):
        barcode = list(self.data.keys())
        barcode_to_id = {barcode: i + 1 for i, barcode in enumerate(barcode)}
        feature = list(set([feature for barcode in self.data for feature in self.data[barcode]]))
        feature_to_id = {feature: i + 1 for i, feature in enumerate(feature)}
        
        with open(self.barcode_file, 'w') as f:
            for bc in barcode:
                f.write(bc + '\n')
        
        with open(self.feature_file, 'w') as f:
            for ft in feature:
                f.write(ft + '\n')
        
        with open(self.matrix_file, 'w') as f:
            content = []
            for bc in barcode:
                for ft in feature:
                    if ft in self.data[bc]:
                        content.append([feature_to_id[ft], barcode_to_id[bc], self.data[bc][ft]])
            
            f.write('%%MatrixMarket matrix coordinate integer general\n')
            f.write('%metadata_json: {"software_version": "cellranger-8.0.1", "format_version": 2}\n')
            f.write(f'{len(feature)} {len(barcode)} {len(content)}\n')
            for line in content:
                f.write(' '.join([str(i) for i in line]) + '\n')


if __name__ == '__main__':
    data_paths = {
        'species1': {
            'barcode': r'D:\Bioinformatics\数据整合\最终版数据\下游分析\SP-data\barcodes.tsv',
            'feature': r'D:\Bioinformatics\数据整合\最终版数据\下游分析\SP-data\features.tsv',
            'matrix': r'D:\Bioinformatics\数据整合\最终版数据\下游分析\SP-data\matrix.mtx'
        },
        'species2': {
            'barcode': r'D:\Bioinformatics\数据整合\最终版数据\下游分析\PD-data\barcodes.tsv',
            'feature': r'D:\Bioinformatics\数据整合\最终版数据\下游分析\PD-data\features.tsv',
            'matrix': r'D:\Bioinformatics\数据整合\最终版数据\下游分析\PD-data\matrix.mtx'
        },
        'species3': {
            'barcode': r'D:\Bioinformatics\数据整合\最终版数据\下游分析\PV-data\barcodes.tsv',
            'feature': r'D:\Bioinformatics\数据整合\最终版数据\下游分析\PV-data\features.tsv',
            'matrix': r'D:\Bioinformatics\数据整合\最终版数据\下游分析\PV-data\matrix.mtx'
        },
        'species4': {
            'barcode': r'D:\nextcloud\pd论文\data\merge-data\FL-data\barcodes.tsv',
            'feature': r'D:\nextcloud\pd论文\data\merge-data\FL-data\features.tsv',
            'matrix': r'D:\nextcloud\pd论文\data\merge-data\FL-data\matrix.mtx'
        },
        'species5': {
            'barcode': r'D:\nextcloud\pd论文\data\merge-data\SD-data\barcodes.tsv',
            'feature': r'D:\nextcloud\pd论文\data\merge-data\SD-data\features.tsv',
            'matrix': r'D:\nextcloud\pd论文\data\merge-data\SD-data\matrix.mtx'
        },
        'species6':{
            'barcode': r'D:\nextcloud\pd论文\data\merge-data\SP-NCBI-data\barcodes.tsv',
            'feature': r'D:\nextcloud\pd论文\data\merge-data\SP-NCBI-data\features.tsv',
            'matrix': r'D:\nextcloud\pd论文\data\merge-data\SP-NCBI-data\matrix.mtx'
        }
    }
    
    merge_data_list = []
    
    for species, paths in data_paths.items():
        try:
            print(f"Processing {species}...")
            barcode = Barcode(paths['barcode'])
            feature = Feature(paths['feature'])
            matrix = Matrix(paths['matrix'])
            
            data = MergeData(barcode, feature, matrix)
            merge_data_list.append(data)
            print(f"Successfully processed {species}")
        except Exception as e:
            print(f"Error processing {species}: {str(e)}")
            continue
    
    if not merge_data_list:
        print("No data was successfully processed. Exiting...")
        exit()
    
    final_data = merge_data_list[0]
    for data in merge_data_list[1:]:
        final_data.merge_data(data)
    
    out_file_path = r"D:\nextcloud\pd论文\result\merge_data"
    os.makedirs(out_file_path, exist_ok=True)
    save_data = SaveData(final_data.data, out_file_path)
    
    print('Merging completed successfully')
