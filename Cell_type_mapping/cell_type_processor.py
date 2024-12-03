from read_marker import ReadMarker
from read_ortho_group import ReadOrthoGroup
from mapping import Mapping

class CellTypeProcessor:
    def __init__(self):
        self.results = {}
        
    def process_species_pairs(self, data_pairs):
        """处理所有物种对
        
        Args:
            data_pairs (dict): 包含物种对文件路径的字典
            
        Returns:
            dict: 处理结果字典，格式为 {species_pair: mapping_results}
        """
        results = {}
        for name, files in data_pairs.items():
            print(f"\n处理 {name} 物种对...")
            mapping = self.process_single_pair(
                files['ortho'],
                files['marker1'],
                files['marker2']
            )
            # 直接获取映射结果而不是返回Mapping对象
            results[name] = mapping.get_mapping_results()
        return results

    @staticmethod
    def process_single_pair(ortho_file, marker1_file, marker2_file):
        """处理单个物种对
        
        Args:
            ortho_file (str): 正交群文件路径
            marker1_file (str): 第一个物种的marker文件路径
            marker2_file (str): 第二个物种的marker文件路径
            
        Returns:
            Mapping: 映射对象
        """
        ortho_group = ReadOrthoGroup(ortho_file, sp1_col=2, sp2_col=1)
        marker1 = ReadMarker(marker1_file)
        marker2 = ReadMarker(marker2_file)

        mapping = Mapping(ortho_group, marker1, marker2)
        mapping.test_overlap()
        mapping.mapping(min_overlap=3, min_confidence=0.05)
        
        return mapping