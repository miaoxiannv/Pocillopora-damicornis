"""
File GSEA.py
Description: 
   This module provides a GSEA (Gene Set Enrichment Analysis) class for performing gene set enrichment analysis.
   It supports analysis using various annotation databases such as GO, KEGG, KOG, and Pfam.
Author: Shengyao Zhang
ate: 2024-12-19
Usage:
   1. Create a GSEA object by providing the gene list file and annotation files.
   2. Run the analysis using the `run_analysis()` method.
   3. Access the analysis results through the corresponding result attributes (e.g., `go_res`, `kegg_res`, etc.).
   4. Optionally, save the analysis results to files using the `save_result()` method.
Example:
   gene_list_file = "path/to/gene_list.txt"
   anno_files = {
       "go": "path/to/go_annotation.xls",
       "kegg": "path/to/kegg_annotation.xls",
       "kog": "path/to/kog_annotation.xls",
       "pfam": "path/to/pfam_annotation.xls"
   }
   
   gsea = GSEA(gene_list_file, **anno_files)
   go_res, go_bp_res, go_cc_res, go_mf_res, kegg_res, kog_res, pfam_res = gsea.run_analysis()
   
   gsea.save_result("path/to/output_directory", format='tsv')
Dependencies:
   - data.gene_list_obj
   - data.gene_set_obj_go_trans
   - data.gene_set_obj_kegg
   - data.gene_set_obj_kog_trans
   - data.gene_set_obj_pfam_trans
   - algorithms.hypergeom
   - os
   - pandas
"""
from data.gene_list_obj import GeneList_Obj
from data import gene_set_obj_go_trans
from data import gene_set_obj_kegg
from data import gene_set_obj_kog_trans
from data import gene_set_obj_pfam_trans
from algorithms import hypergeom
import os
import pandas as pd

class GSEA:
    def __init__(self, gene_list_file, **anno_files):
        # Initialize file paths
        self.gene_list_file = gene_list_file
        self.anno_files = anno_files
        

        # Initialize result variables
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
        
        # Import gene list
        self.gene_list = GeneList_Obj(self.gene_list_file)
        
        # Import annotation files
        self._load_annotations()
    
    def _load_annotations(self):
        """Load all annotation files"""

        
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
        """Run GSEA analysis"""
        # Run hypergeometric test
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
        Save analysis results
        Args:
            out_dir: Output directory path
            format: Output format, supports 'xlsx' or 'tsv'
        """
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        # Define save function
        def save_df(df, filename):
            if df is not None and not df.empty:
                file_path = os.path.join(out_dir, filename)
                if format.lower() == 'xlsx':
                    df.to_excel(f"{file_path}.xlsx", index=False)
                elif format.lower() == 'tsv':
                    df.to_csv(f"{file_path}.tsv", sep='\t', index=False)
                else:
                    raise ValueError(f"Unsupported format: {format}")
        
        # Save various results
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
    """Test main function for GSEA analysis"""
    # Test file paths
    gene_list_file = r"D:\nextcloud\pd论文\data\Cluster-pos-nag-genlist\Cluster0_markers_negative.txt"
    anno_files = {
        "go": r"D:\nextcloud\pd论文\data\test\gsea\GO.anno.xls",
        "kegg": r"D:\nextcloud\pd论文\data\test\gsea\P2.KEGG.filter.m8.anno.xls",
        "kog": r"D:\nextcloud\pd论文\data\test\gsea\P2.KOG.filter.m8.anno.xls",
        "pfam": r"D:\nextcloud\pd论文\data\test\gsea\P2.pfam.anno.xls"
    }

    try:
        # Create GSEA object
        gsea = GSEA(gene_list_file, **anno_files)
        
        # Run analysis
        go_res, go_bp_res, go_cc_res, go_mf_res, kegg_res, kog_res, pfam_res = gsea.run_analysis()
        
        # Print results
        print("\nGO Enrichment Analysis Results:")
        if go_res is not None:
            print(f"Total GO terms: {len(go_res)}")
            print(go_res.head())
            
            print("\nGO Subcategory Analysis Results:")
            print(f"BP terms: {len(go_bp_res)}")
            print(go_bp_res.head())
            print(f"CC terms: {len(go_cc_res)}")
            print(go_cc_res.head())
            print(f"MF terms: {len(go_mf_res)}")
            print(go_mf_res.head())

        
        print("\nKEGG Enrichment Analysis Results:")
        if kegg_res is not None:
            print(f"Enriched KEGG terms: {len(kegg_res)}")
            print(kegg_res.head())
            
        print("\nKOG Enrichment Analysis Results:")
        if kog_res is not None:
            print(f"Enriched KOG terms: {len(kog_res)}")
            print(kog_res.head())
            
        print("\nPfam Enrichment Analysis Results:")
        if pfam_res is not None:
            print(f"Enriched Pfam terms: {len(pfam_res)}")
            print(pfam_res.head())
            
    except Exception as e:
        print(f"Error during execution: {str(e)}")

    # Export result
    gsea.save_result(r"D:\nextcloud\pd论文\result\GSEA\Cluster0\negative", format='tsv')




