from collections import defaultdict
import os


class Mapping:
    def __init__(self, ortho_group, marker_sp1, marker_sp2):
        self.ortho_group = ortho_group
        self.gene_mapping_dict_sp1 = self.ortho_group.get_gene_mapping_dict_sp1()
        self.gene_mapping_dict_sp2 = self.ortho_group.get_gene_mapping_dict_sp2()
        self.marker_sp1 = marker_sp1
        self.marker_sp2 = marker_sp2
        self.marker_dict_sp1 = self.marker_sp1.get_marker()
        self.marker_dict_sp2 = self.marker_sp2.get_marker()

    def test_overlap(self):
        """测试基因重叠情况并输出统计信息"""
        sp1_marker_genes = self._collect_genes(self.marker_dict_sp1)
        sp2_marker_genes = self._collect_genes(self.marker_dict_sp2)
        sp1_ortho_genes = self.ortho_group.sp1_all_genes
        sp2_ortho_genes = self.ortho_group.sp2_all_genes

        print("基因统计:")
        print(f"sp1 marker基因数量: {len(sp1_marker_genes)}")
        print(f"sp2 marker基因数量: {len(sp2_marker_genes)}")
        print(f"sp1 直系同源基因数量: {len(sp1_ortho_genes)}")
        print(f"sp2 直系同源基因数量: {len(sp2_ortho_genes)}")
        print(f"sp1 marker基因与直系同源基因重叠数量: {len(sp1_marker_genes.intersection(sp1_ortho_genes))}")
        print(f"sp2 marker基因与直系同源基因重叠数量: {len(sp2_marker_genes.intersection(sp2_ortho_genes))}")

    def mapping(self, min_overlap=3, min_confidence=0.01):
        """
        执行细胞类型映射

        参数:
        min_overlap: 最小重叠基因数量
        min_confidence: 最小置信度阈值 (0-1之间)
        """
        mapping_results = defaultdict(list)

        for cell_type_sp2, genes_sp2 in self.marker_dict_sp2.items():
            cell_mappings = self._map_cell_type(cell_type_sp2, genes_sp2, min_overlap, min_confidence)
            mapping_results[cell_type_sp2] = cell_mappings

        self._print_mapping_results(mapping_results)

    def _map_cell_type(self, cell_type_sp2, genes_sp2, min_overlap, min_confidence):
        """映射单个细胞类型"""
        cell_mappings = []

        for cell_type_sp1, genes_sp1 in self.marker_dict_sp1.items():
            transformed_genes_sp1 = self.get_transformed_genes(genes_sp1, self.gene_mapping_dict_sp1)
            overlap_genes = transformed_genes_sp1.intersection(genes_sp2)

            if len(overlap_genes) >= min_overlap:
                confidence = len(overlap_genes) / len(genes_sp2)
                if confidence >= min_confidence:
                    cell_mappings.append({
                        'sp1_cell_type': cell_type_sp1,
                        'overlap_count': len(overlap_genes),
                        'confidence': confidence,
                        'overlap_genes': overlap_genes
                    })

        cell_mappings.sort(key=lambda x: x['confidence'], reverse=True)
        return cell_mappings

    def _print_mapping_results(self, mapping_results):
        """打印映射结果"""
        for cell_type_sp2, mappings in mapping_results.items():
            print(f"\n细胞类型 {cell_type_sp2} 的映射结果:")
            if mappings:
                for mapping in mappings:
                    print(f"-> {mapping['sp1_cell_type']}")
                    print(f"   重叠基因数: {mapping['overlap_count']}")
                    print(f"   置信度: {mapping['confidence']:.2%}")
                    print(f"   重叠基因: {', '.join(sorted(mapping['overlap_genes']))}")
                    print()
            else:
                print("未找到满足条件的映射关系")

    def get_transformed_genes(self, genes, gene_mapping_dict):
        """转换基因名称"""
        gene_set = set()
        for gene in genes:
            if gene in gene_mapping_dict:
                gene_set.update(gene_mapping_dict[gene])
        return gene_set


    @staticmethod
    def _collect_genes(marker_dict):
        """收集所有细胞类型的基因"""
        all_genes = set()
        for genes in marker_dict.values():
            all_genes.update(genes)
        return all_genes

    def get_statistics(self):
        """获取统计信息"""
        sp1_marker_genes = self._collect_genes(self.marker_dict_sp1)
        sp2_marker_genes = self._collect_genes(self.marker_dict_sp2)
        sp1_ortho_genes = self.ortho_group.sp1_all_genes
        sp2_ortho_genes = self.ortho_group.sp2_all_genes
        
        sp1_overlap = len(sp1_marker_genes.intersection(sp1_ortho_genes))
        sp2_overlap = len(sp2_marker_genes.intersection(sp2_ortho_genes))
        
        return {
            'sp1_marker': len(sp1_marker_genes),
            'sp2_marker': len(sp2_marker_genes),
            'sp1_ortho': len(sp1_ortho_genes),
            'sp2_ortho': len(sp2_ortho_genes),
            'sp1_overlap': sp1_overlap,
            'sp2_overlap': sp2_overlap,
            'sp1_overlap_rate': sp1_overlap / len(sp1_marker_genes),
            'sp2_overlap_rate': sp2_overlap / len(sp2_marker_genes)
        }

    def get_mapping_results(self):
        """获取映射结果"""
        mapping_results = {}
        for cell_type_sp2, genes_sp2 in self.marker_dict_sp2.items():
            mappings = self._map_cell_type(cell_type_sp2, genes_sp2, 3, 0.01)
            if mappings:
                mapping_results[cell_type_sp2] = mappings
        return mapping_results
