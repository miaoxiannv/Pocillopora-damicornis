"""
Gene list consistency checker.
This script compares two gene lists and identifies overlapping genes between them,
helping ensure consistency between different gene sets.

Usage:
    python ensure_consistency.py
    This will compare genes between two specified files and report overlapping genes.

Author: Shengyao Zhang
Date: 2024-12-19
Version: 1.0
"""

def read_genes_from_file(file_path):
    """Read gene names from file and return a set"""
    with open(file_path, 'r') as f:
        genes = {line.strip() for line in f if line.strip()}
    return genes

def find_overlapping_genes(file1, file2):
    """Find overlapping genes between two files"""
    genes1 = read_genes_from_file(file1)
    genes2 = read_genes_from_file(file2)
    
    # Find overlapping genes
    overlapping_genes = genes1.intersection(genes2)
    
    print(f"Number of genes in file 1: {len(genes1)}")
    print(f"Number of genes in file 2: {len(genes2)}")
    print(f"Number of overlapping genes: {len(overlapping_genes)}")
    
    return overlapping_genes

if __name__ == "__main__":
    file1 = r"D:\nextcloud\pd论文\data\test\Cell_type_prediction\PD-SP-cell\PD-SP-cell_Orthogroups_genes.txt"  # Replace with path to first file
    file2 = r"D:\nextcloud\pd论文\data\test\Cell_type_prediction\PD-SP-cell\feature_names.tsv"  # Replace with path to second file
    
    # Find overlapping genes
    overlapping_genes = find_overlapping_genes(file1, file2)
    
    # Output overlapping genes
    if overlapping_genes:
        print("Number of overlapping genes:")
        print(len(overlapping_genes))
