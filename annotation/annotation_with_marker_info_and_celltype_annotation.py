"""
Annotate marker genes with cell type information.

This script reads a marker gene file and an annotation file, merges the information,
and exports the annotated marker genes to an output file.

Usage:
    python annotation_with_marker_info_and_celltype_annotation.py

Author: Shengyao Zhang
Date: 2024-12-19
Version: 1.0
"""

import os
import pandas as pd

class AnnotationCellType:
    """
    Class for annotating marker genes with cell type information.
    """
    def __init__(self, marker_file, annotation_file):
        """
        Initialize the AnnotationCellType with marker file and annotation file paths.
        
        Args:
            marker_file (str): Path to the marker gene file.
            annotation_file (str): Path to the annotation file.
        """
        self.marker_file = marker_file
        self.annotation_file = annotation_file
        self.marker_df = None
        self.annotation_df = None
        self.merged_df = None
    
    def read_marker_file(self):
        """
        Read the marker gene file, including p_val and other statistical information.
        
        Returns:
            pd.DataFrame: The marker gene DataFrame.
        """
        # Read the marker file, including p_val and other statistical information
        self.marker_df = pd.read_csv(self.marker_file, sep='\t', header=0)
        
        # Ensure required statistical columns exist
        required_columns = ['Gene_id', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj']
        for col in required_columns:
            if col not in self.marker_df.columns:
                print(f"Warning: Missing column {col}")
                self.marker_df[col] = ''
        
        # Normalize gene IDs
        self.marker_df['Gene_id'] = self.marker_df['Gene_id'].str.strip().str.lower()
        return self.marker_df
    
    def read_annotation_file(self):
        """
        Read the annotation file, without using the default header.
        
        Returns:
            pd.DataFrame: The annotation DataFrame.
        """
        # Read the annotation file, without using the default header
        self.annotation_df = pd.read_csv(self.annotation_file, sep='\t', header=None)
        
        # Set column names based on file content
        # First column is Gene_id, second column is Domain, remaining columns are merged into Description
        columns = ['Gene_id', 'Domain']
        for i in range(2, len(self.annotation_df.columns)):
            columns.append(f'Description{i-1}')
        
        self.annotation_df.columns = columns
        
        # Normalize gene IDs
        self.annotation_df['Gene_id'] = self.annotation_df['Gene_id'].str.strip().str.lower()
        
        # Clean the Domain column (remove trailing slash)
        self.annotation_df['Domain'] = self.annotation_df['Domain'].str.rstrip('/')
        
        return self.annotation_df
    
    def merge_data(self):
        """
        Merge marker gene information and annotation information.
        
        Returns:
            pd.DataFrame: The merged DataFrame.
        """
        # Merge marker information and annotation information
        self.merged_df = pd.merge(self.marker_df, 
                                self.annotation_df,
                                on='Gene_id',
                                how='left')
        
        # Statistics of the merge result
        total_genes = len(self.marker_df)
        annotated_genes = self.merged_df['Domain'].notna().sum()
        print(f"Total genes: {total_genes}")
        print(f"Annotated genes: {annotated_genes}")
        print(f"Unannotated genes: {total_genes - annotated_genes}")
        
        return self.merged_df
    
    def export(self, out_file):
        """
        Export the annotated marker genes to an output file.
        
        Args:
            out_file (str): Path to the output file.
        """
        # Ensure the output directory exists
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        
        # Sort by p-value and export
        sorted_df = self.merged_df.sort_values(by='p_val', ascending=True)
        
        # Export the file
        if out_file.endswith('.xlsx'):
            sorted_df.to_excel(out_file, index=False)
        else:
            sorted_df.to_csv(out_file, sep='\t', index=False)
        
        print(f"Annotation file saved to: {out_file}")

if __name__ == "__main__":
    # Input and output file paths
    marker_file = r"D:\nextcloud\pd论文\data\cluster-marker\cluster0_markers.tsv"
    annotation_file = r"D:\nextcloud\pd论文\data\g2mysijfp52-3\Spis_gene_annotation.tsv"  # Gene annotation file
    output_file = r"D:\nextcloud\pd论文\result\annotation1\cluster0_markers_annotated.tsv"
    
    # Process the annotation
    anno = AnnotationCellType(marker_file, annotation_file)
    anno.read_marker_file()
    anno.read_annotation_file()
    anno.merge_data()
    anno.export(output_file)























