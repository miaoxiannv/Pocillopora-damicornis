'''
1. read gene list;
2. read go, kegg, kog and pfam annotation xlsx file;
3. convert gene list to go, kegg, kog and pfam;
4. run hypergeometric test;
5. export result;
6. plot result;
'''
#import gene list object
from data.gene_list_obj import GeneList_Obj
#import gene set object
from data import gene_set_obj_go_trans
from data import gene_set_obj_kegg
from data import gene_set_obj_kog_trans
from data import gene_set_obj_pfam_trans

#import hypergeometric test
from algorithms import hypergeom

#import plot
from evaluation import gsea_plot_no_type
from evaluation import gsea_plot_with_type


import os
import pandas as pd

class GSEA:
    def __init__(self, gene_list_file, **anno_files):
        # 初始化文件路径
        self.gene_list_file = gene_list_file
        self.anno_files = anno_files
        

        # 初始化结果变量
        self.go_gmt = None
        self.go_sub_objs = None
        self.kegg_gmt = None
        self.kog_gmt = None
        self.pfam_gmt = None

        self.go_res = None
        self.go_bp_res = None
        self.go_cc_res = None
        self.go_mf_res = None
        self.kegg_res = None
        self.kog_res = None
        self.pfam_res = None
        
        # 导入基因列表
        self.gene_list = GeneList_Obj(self.gene_list_file)
        
        # 导入注释文件
        self._load_annotations()
    
    def _load_annotations(self):
        """加载所有注释文件"""

        
        for gene_set, file_path in self.anno_files.items():
            if gene_set == 'go':
                self.go_gmt = gene_set_obj_go_trans.Gmt_stat(file_path)
                self.go_sub_objs = self.go_gmt.split_by_term_class()
            elif gene_set == 'kegg':
                self.kegg_gmt = gene_set_obj_kegg.Gmt_stat_gene(file_path)
            elif gene_set == 'kog':
                self.kog_gmt = gene_set_obj_kog_trans.Gmt_stat(file_path)
            elif gene_set == 'pfam':
                self.pfam_gmt = gene_set_obj_pfam_trans.Gmt_stat(file_path)
                            
    def run_analysis(self):
        """运行GSEA分析"""
        # 运行超几何检验
        if self.go_gmt:
            self.go_res = hypergeom.calcu_hypergeom(self.gene_list, self.go_gmt)
            for class_type, sub_obj in self.go_sub_objs.items():
                if class_type.upper() == 'BP':
                    self.go_bp_res = hypergeom.calcu_hypergeom(self.gene_list, sub_obj)
                elif class_type.upper() == 'CC':
                    self.go_cc_res = hypergeom.calcu_hypergeom(self.gene_list, sub_obj)
                elif class_type.upper() == 'MF':
                    self.go_mf_res = hypergeom.calcu_hypergeom(self.gene_list, sub_obj)
        
        if self.kegg_gmt:
            self.kegg_res = hypergeom.calcu_hypergeom(self.gene_list, self.kegg_gmt)
        
        if self.kog_gmt:
            self.kog_res = hypergeom.calcu_hypergeom(self.gene_list, self.kog_gmt)
        
        if self.pfam_gmt:
            self.pfam_res = hypergeom.calcu_hypergeom(self.gene_list, self.pfam_gmt)
        
        return self.go_res, self.go_bp_res, self.go_cc_res, self.go_mf_res, self.kegg_res, self.kog_res, self.pfam_res
    
    #save all result
    def save_result(self, out_dir, format='xlsx'):
        """
        保存分析结果
        Args:
            out_dir: 输出目录路径
            format: 输出格式，支持 'xlsx' 或 'tsv'
        """
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        # 定义保存函数
        def save_df(df, filename):
            if df is not None and not df.empty:
                file_path = os.path.join(out_dir, filename)
                if format.lower() == 'xlsx':
                    df.to_excel(f"{file_path}.xlsx", index=False)
                elif format.lower() == 'tsv':
                    df.to_csv(f"{file_path}.tsv", sep='\t', index=False)
                else:
                    raise ValueError(f"不支持的格式: {format}")
        
        # 保存各类结果
        result_files = {
            'go_res': self.go_res,
            'go_bp_res': self.go_bp_res,
            'go_cc_res': self.go_cc_res,
            'go_mf_res': self.go_mf_res,
            'kegg_res': self.kegg_res,
            'kog_res': self.kog_res,
            'pfam_res': self.pfam_res
        }
        
        for name, df in result_files.items():
            save_df(df, name)


if __name__ == "__main__":
    """测试GSEA分析的主函数"""
    # 测试文件路径
    gene_list_file = r"D:\nextcloud\pd论文\data\Cluster-pos-nag-genlist\cluster0_markers_negative.txt"
    anno_files = {
        "go": r"D:\nextcloud\pd论文\data\test\gsea\GO.anno.xls",
        "kegg": r"D:\nextcloud\pd论文\data\test\gsea\P2.KEGG.filter.m8.anno.xls",
        "kog": r"D:\nextcloud\pd论文\data\test\gsea\P2.KOG.filter.m8.anno.xls",
        "pfam": r"D:\nextcloud\pd论文\data\test\gsea\P2.pfam.anno.xls"
    }

    try:
        # 创建GSEA对象
        gsea = GSEA(gene_list_file, **anno_files)
        
        # 运行分析
        go_res, go_bp_res, go_cc_res, go_mf_res, kegg_res, kog_res, pfam_res = gsea.run_analysis()
        
        # 打印结果
        print("\nGO富集分析结果:")
        if go_res is not None:
            print(f"GO总条目数: {len(go_res)}")
            print(go_res.head())
            
            print("\nGO子类别分析结果:")
            print(f"BP条目数: {len(go_bp_res)}")
            print(go_bp_res.head())
            print(f"CC条目数: {len(go_cc_res)}")
            print(go_cc_res.head())
            print(f"MF条目数: {len(go_mf_res)}")
            print(go_mf_res.head())

        
        print("\nKEGG富集分析结果:")
        if kegg_res is not None:
            print(f"KEGG富集条目数: {len(kegg_res)}")
            print(kegg_res.head())
            
        print("\nKOG富集分析结果:")
        if kog_res is not None:
            print(f"KOG富集条目数: {len(kog_res)}")
            print(kog_res.head())
            
        print("\nPfam富集分析结果:")
        if pfam_res is not None:
            print(f"Pfam富集条目数: {len(pfam_res)}")
            print(pfam_res.head())
            
    except Exception as e:
        print(f"运行过程中出错: {str(e)}")

    # export result
    gsea.save_result(r"D:\nextcloud\pd论文\result\CollateCode\GSEA\Cluster0\negative", format='tsv')




