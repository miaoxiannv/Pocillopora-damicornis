"""
Annotte marker genes with additional information from various databases.
This script reads a marker gene file and multiple annotation files (GO, KEGG, Pfam, KOG),
erges the information, and exports the annotated marker genes to output files.
Usage:
   python annotation_with_marker_info.py
Author: Shengyao Zhang
ate: 2024-12-19
ersion: 1.0
"""
import os
import pandas as pd
from abc import ABC, abstractmethod
class GeneList:
   """
   Class for reading and processing a gene list file.
   """
   def __init__(self, gene_list_file):
       """
       Initialize the GeneList with the gene list file path.
       
       Args:
           gene_list_file (str): Path to the gene list file.
       """
       self.gene_list_file = gene_list_file
       self.gene_list = None
       self.read_gene_list()
   
   def read_gene_list(self):
       """
       Read the gene list file and process the gene names.
       """
       gene_df = pd.read_csv(self.gene_list_file, header=None, index_col=False)
       gene_df.columns = ['gene']
       gene_df['gene'] = gene_df['gene'].str.strip().str.lower()
       self.gene_list = gene_df['gene'].tolist()
   
   def get_gene_list(self):
       """
       Get the processed gene list.
       
       Returns:
           list: The gene list.
       """
       return self.gene_list
class Annotation(ABC):
   """
   Abstract base class for annotation classes.
   """
   def __init__(self, anno_file):
       """
       Initialize the Annotation with the annotation file path.
       
       Args:
           anno_file (str): Path to the annotation file.
       """
       self.anno_file = anno_file
       self.anno_df = None
   
   def read_anno_file(self):
       """
       Read the annotation file and process the data.
       First read the header line to get column names.
       """
       # First try to read the first line to get column names
       with open(self.anno_file, 'r') as f:
           header = f.readline().strip().split('\t')
       
       # Read the entire file content
       data = []
       with open(self.anno_file, 'r') as f:
           next(f)  # Skip header line
           for line in f:
               fields = line.strip().split('\t')
               if len(fields) > len(header):
                   # If number of fields exceeds header count, merge extra fields into last column
                   merged_fields = fields[:len(header)-1] + ['\t'.join(fields[len(header)-1:])]
                   data.append(merged_fields)
               else:
                   data.append(fields)
       
       # Create DataFrame
       self.anno_df = pd.DataFrame(data, columns=header)
       return self.anno_df
   
   @abstractmethod
   def set_index(self):
       """
       Abstract method for setting the index of the annotation DataFrame.
       """
       pass
   
   def search(self, gene_list):
       """
       Search for the annotations of the given gene list.
       
       Args:
           gene_list (list): The list of genes to search for.
       
       Returns:
           pd.DataFrame: The search results DataFrame.
       """
       gene_list = [gene.strip().lower() for gene in gene_list]
       
       res = self.anno_df.reindex(gene_list)
       found_genes = res.dropna(how='all').shape[0]
       res.reset_index(inplace=True)
       
       not_found_genes = res.shape[0] - found_genes
       
       print(f"Search Results:")
       print(f"- Total input genes: {len(gene_list)}")
       print(f"- Genes with annotation: {found_genes}")
       print(f"- Genes without annotation: {not_found_genes}")
       
       return res
   
   def export(self, res, out_file):
       """
       Export the search results to an output file.
       
       Args:
           res (pd.DataFrame): The search results DataFrame.
           out_file (str): Path to the output file.
       """
       # Create base filename (without extension)
       base_file = os.path.splitext(out_file)[0]
       
       # Export original results
       if out_file.endswith('.xlsx'):
           res.to_excel(out_file, index=False)
           # Get second column name
           second_col = res.columns[1]
           # Create sorted results
           sorted_res = res.sort_values(by=second_col, ascending=False)
           # Export sorted results
           sorted_res.to_excel(f"{base_file}.sorted.xlsx", index=False)
       elif out_file.endswith('.tsv'):
           res.to_csv(out_file, sep='\t', index=False)
           # Get second column name
           second_col = res.columns[1]
           # Create sorted results
           sorted_res = res.sort_values(by=second_col, ascending=False)
           # Export sorted results
           sorted_res.to_csv(f"{base_file}.sorted.tsv", sep='\t', index=False)
       else:
           raise ValueError("Unsupported file format. Please use .xlsx or .tsv.")
class AnnotationGo(Annotation):
   """
   Class for GO annotation.
   """
   def set_index(self):
       """
       Set the index of the GO annotation DataFrame to 'Gene_id'.
       """
       self.anno_df['Gene_id'] = self.anno_df['Gene_id'].str.strip().str.lower()
       self.anno_df.set_index('Gene_id', inplace=True)
class AnnotationKegg(Annotation):
   """
   Class for KEGG annotation.
   """
   def set_index(self):
       """
       Set the index of the KEGG annotation DataFrame to 'Query_id'.
       """
       self.anno_df['Query_id'] = self.anno_df['Query_id'].str.strip().str.lower()
       
       # Check for duplicate records
       duplicates = self.anno_df['Query_id'].duplicated().sum()
       print(f"Found {duplicates} duplicate gene IDs")
       
       # Group and merge duplicate records
       # Merge all columns except Query_id
       merged_df = self.anno_df.groupby('Query_id').agg({
           'Subject_id': lambda x: ' || '.join(x.dropna().unique()),
           'KO_ID': lambda x: ' || '.join(x.dropna().unique()),
           'KO_NAME': lambda x: ' || '.join(x.dropna().unique()),
           'KO_DEFINITION': lambda x: ' || '.join(x.dropna().unique()),
           'KO_EC': lambda x: ' || '.join(x.dropna().unique()),
           'KO_PATHWAY': lambda x: ' || '.join(x.dropna().unique())
       }).reset_index()
       
       # Set index
       self.anno_df = merged_df.set_index('Query_id')
class AnnotationPfam(Annotation):
   """
   Class for Pfam annotation.
   """
   def set_index(self):
       """
       Set the index of the Pfam annotation DataFrame to 'Gene_id'.
       """
       self.anno_df['Gene_id'] = self.anno_df['Gene_id'].str.strip().str.lower()
       self.anno_df.set_index('Gene_id', inplace=True)
class AnnotationKog(Annotation):
   """
   Class for KOG annotation.
   """
   def set_index(self):
       """
       Set the index of the KOG annotation DataFrame to 'Gene_id'.
       """
       self.anno_df['Gene_id'] = self.anno_df['Gene_id'].str.strip().str.lower()
       self.anno_df.set_index('Gene_id', inplace=True)
class AnnotationMarker(Annotation):
   """
   Class for marker gene annotation.
   """
   def read_anno_file(self):
       """
       Read the marker file with pandas directly
       """
       # Read tsv file directly using pandas
       self.anno_df = pd.read_csv(self.anno_file, sep='\t', header=0)
       
       # Print column names for debugging
       print("Marker file columns:", self.anno_df.columns.tolist())
       
       # Rename first column to 'gene' if unnamed
       if self.anno_df.columns[0].startswith('Unnamed'):
           self.anno_df.rename(columns={self.anno_df.columns[0]: 'gene'}, inplace=True)
       
       # Replace hyphens with underscores in gene column
       self.anno_df[self.anno_df.columns[0]] = self.anno_df[self.anno_df.columns[0]].str.replace('-', '_')
       
       return self.anno_df
   
   def set_index(self):
       """
       Set the index of the marker gene annotation DataFrame.
       """
       # Get first column as gene name column
       gene_col = self.anno_df.columns[0]
       
       # Clean and convert gene names to lowercase
       self.anno_df[gene_col] = self.anno_df[gene_col].str.strip().str.lower()
       
       # Check for duplicate records
       duplicates = self.anno_df[gene_col].duplicated().sum()
       if duplicates > 0:
           print(f"Found {duplicates} duplicate gene IDs")
           # Group and merge duplicates, keeping the record with smallest p-value
           self.anno_df = self.anno_df.sort_values('p_val').groupby(gene_col).first().reset_index()
       
       # Set index
       self.anno_df.set_index(gene_col, inplace=True)
if __name__ == "__main__":
   # Input files
   marker_file = r"D:\nextcloud\pd论文\data\cluster-marker\cluster9_markers.tsv"
   go_anno_file = r"D:\nextcloud\pd论文\data\test\annotation\GO.anno.xls"
   kegg_anno_file = r"D:\nextcloud\pd论文\data\test\annotation\P2.KEGG.filter.m8.anno.xls"
   pfam_anno_file = r"D:\nextcloud\pd论文\data\test\annotation\P2.pfam.anno.xls"
   kog_anno_file = r"D:\nextcloud\pd论文\data\test\annotation\P2.KOG.filter.m8.anno.xls"
   
   # Output directory
   out_dir = r"D:\nextcloud\pd论文\result\annotation\cluster9"
   if not os.path.exists(out_dir):
       os.makedirs(out_dir)
   
   # First read and process marker file
   print("Testing marker gene annotation...")
   marker_anno = AnnotationMarker(marker_file)
   marker_anno.read_anno_file()
   marker_anno.set_index()
   print(f"Number of marker genes: {len(marker_anno.anno_df)}")
   print(marker_anno.anno_df.head())
   print("\n" + "="*50 + "\n")
   
   # Use gene list from marker file
   gene_list = marker_anno.anno_df.index.tolist()
   
   # Test GO annotation
   print("Testing GO annotation...")
   go_anno = AnnotationGo(go_anno_file)
   go_anno.read_anno_file()
   go_anno.set_index()
   go_results = go_anno.search(gene_list)
   go_results = pd.merge(go_results, 
                        marker_anno.anno_df.reset_index(), 
                        left_on='Gene_id',
                        right_on='gene', 
                        how='left')
   
   print(f"Number of GO annotation results: {len(go_results)}")
   print(go_results.head())
   go_anno.export(go_results, os.path.join(out_dir, "go_results.tsv"))
   print("\n" + "="*50 + "\n")
    # Test KEGG annotation
   print("Testing KEGG annotation...")
   kegg_anno = AnnotationKegg(kegg_anno_file)
   kegg_anno.read_anno_file()
   kegg_anno.set_index()
   kegg_results = kegg_anno.search(gene_list)
   kegg_results = pd.merge(kegg_results, 
                          marker_anno.anno_df.reset_index(), 
                          left_on='Query_id',
                          right_on='gene', 
                          how='left')
   print(f"Number of KEGG annotation results: {len(kegg_results)}")
   print(kegg_results.head())
   kegg_anno.export(kegg_results, os.path.join(out_dir, "kegg_results.tsv"))
   print("\n" + "="*50 + "\n")
    # Test Pfam annotation
   print("Testing Pfam annotation...")
   pfam_anno = AnnotationPfam(pfam_anno_file)
   pfam_anno.read_anno_file()
   pfam_anno.set_index()
   pfam_results = pfam_anno.search(gene_list)
   pfam_results = pd.merge(pfam_results, 
                          marker_anno.anno_df.reset_index(), 
                          left_on='Gene_id',
                          right_on='gene', 
                          how='left')
   print(f"Number of Pfam annotation results: {len(pfam_results)}")
   print(pfam_results.head())
   pfam_anno.export(pfam_results, os.path.join(out_dir, "pfam_results.tsv"))
   print("\n" + "="*50 + "\n")
    # Test KOG annotation
   print("Testing KOG annotation...")
   kog_anno = AnnotationKog(kog_anno_file)
   kog_anno.read_anno_file()
   kog_anno.set_index()
   kog_results = kog_anno.search(gene_list)
   kog_results = pd.merge(kog_results, 
                         marker_anno.anno_df.reset_index(), 
                         left_on='Gene_id',
                         right_on='gene', 
                         how='left')
   print(f"Number of KOG annotation results: {len(kog_results)}")
   print(kog_results.head())
   kog_anno.export(kog_results, os.path.join(out_dir, "kog_results.tsv"))
   print("\n" + "="*50 + "\n")