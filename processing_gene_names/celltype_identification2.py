"""
Cell Type Identification Comparator.
This script compares marker genes between two species to identify corresponding cell types.
It processes marker gene files from both species and identifies overlapping genes.

Input files:
    - marker1_files: Marker gene files from first species (LOC111 format)
    - marker2_files: Marker gene files from second species (Spis-XP format)
    - mapping_file: Gene ID to protein ID mapping file

Output:
    - Comparison results for each cell type pair
    - Detailed log files with comparison statistics
    - TSV files containing overlapping genes

Usage:
    python celltype_identification2.py

Author: Shengyao Zhang
Date: 2024-12-19
Version: 1.0
"""

import pandas as pd
from pathlib import Path
from itertools import product
import logging
from datetime import datetime



def read_file_with_encodings(file_path, **kwargs):
    """Generic file reading function that automatically handles encoding issues"""
    encodings = ['utf-8', 'gbk', 'gb18030']
    for encoding in encodings:
        try:
            return pd.read_csv(file_path, encoding=encoding, **kwargs)
        except UnicodeDecodeError:
            continue
    raise UnicodeDecodeError(f"Could not read file with encodings {encodings}: {file_path}")

def process_marker_file1(file_path, gene_to_protein, all_marker1_files, gene_prefix='LOC111', min_log2fc=0, max_pval=0.05):
    """Process the marker file in the first format and retain only the genes that are specifically expressed."""
    # read the current file
    current_markers = read_file_with_encodings(file_path, sep='\t')
    current_proteins = set()
    
    # Read the genes of all other marker1 files
    other_proteins = set()
    for other_file in all_marker1_files:
        if other_file != file_path:  # Skip current file
            other_markers = read_file_with_encodings(other_file, sep='\t')
            for _, row in other_markers.iterrows():
                gene = row.iloc[0]
                avg_log2fc = row['avg_log2FC']
                p_val = row['p_val']
                
                if (gene.startswith(gene_prefix) and 
                    avg_log2fc > min_log2fc and 
                    p_val < max_pval and 
                    gene in gene_to_protein and 
                    gene_to_protein[gene] is not None):
                    protein = gene_to_protein[gene].replace('.', '_')
                    other_proteins.add(protein)
    
    # Only genes that are specifically expressed in the current document are retained.
    for _, row in current_markers.iterrows():
        gene = row.iloc[0]
        avg_log2fc = row['avg_log2FC']
        p_val = row['p_val']
        
        if (gene.startswith(gene_prefix) and 
            avg_log2fc > min_log2fc and 
            p_val < max_pval and 
            gene in gene_to_protein and 
            gene_to_protein[gene] is not None):
            protein = gene_to_protein[gene].replace('.', '_')
            if protein not in other_proteins:  # Add only the proteins that do not appear in other files.
                current_proteins.add(protein)
    
    return list(current_proteins)

def process_marker_file2(file_path, all_marker2_files):
    """Process second type of marker file, keeping only specifically expressed genes"""
    # Read current file
    current_markers = read_file_with_encodings(file_path, sep=" ", quotechar='"')
    current_proteins = set()
    
    # Read genes from all other marker2 files
    other_proteins = set()
    for other_file in all_marker2_files:
        if other_file != file_path:  # Skip current file
            other_markers = read_file_with_encodings(other_file, sep=" ", quotechar='"')
            for gene in other_markers.index:
                # Convert gene to string before using startswith
                gene_str = str(gene)
                if (gene_str.startswith('Spis-XP-') and 
                    other_markers.loc[gene, 'avg_log2FC'] > 0 and 
                    other_markers.loc[gene, 'p_val'] < 0.05):
                    protein = gene_str.replace('Spis-XP-', 'XP_').replace('-', '_')
                    other_proteins.add(protein)
    
    # Keep only genes specifically expressed in current file
    for gene in current_markers.index:
        avg_log2fc = current_markers.loc[gene, 'avg_log2FC']
        p_val = current_markers.loc[gene, 'p_val']
        
        # Convert gene to string before using startswith
        gene_str = str(gene)
        if (gene_str.startswith('Spis-XP-') and 
            avg_log2fc > 0 and 
            p_val < 0.05):
            protein = gene_str.replace('Spis-XP-', 'XP_').replace('-', '_')
            if protein not in other_proteins:  # Only add proteins not in other files
                current_proteins.add(protein)
    
    return list(current_proteins)

def setup_logging(output_dir):
    """Set up logging configuration"""
    log_file = output_dir / f"comparison_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def save_comparison_results(comparison, output_dir):
    """Save comparison results for each pair"""
    # Create result filename
    result_file = output_dir / f"{comparison['marker1_file']}_vs_{comparison['marker2_file']}.tsv"
    
    # Create result dataframe
    result_df = pd.DataFrame({
        'protein_id': sorted(list(comparison['overlapping_proteins'])),
        'source': 'overlapping'
    })
    
    # Add proteins only in marker1
    marker1_only = comparison['marker1_proteins'] - comparison['overlapping_proteins']
    if marker1_only:
        result_df = pd.concat([result_df, pd.DataFrame({
            'protein_id': sorted(list(marker1_only)),
            'source': 'marker1_only'
        })])
    
    # Add proteins only in marker2
    marker2_only = comparison['marker2_proteins'] - comparison['overlapping_proteins']
    if marker2_only:
        result_df = pd.concat([result_df, pd.DataFrame({
            'protein_id': sorted(list(marker2_only)),
            'source': 'marker2_only'
        })])
    
    # Save results
    result_df.to_csv(result_file, sep='\t', index=False)
    return result_file

def compare_marker_pairs(marker1_files, marker2_files, mapping_file, output_dir, 
                        gene_prefix='LOC111', min_log2fc=0, max_pval=0.05):
    """
    Compare overlapping genes between each pair of marker files.
    
    Args:
        marker1_files: List of marker files from first species
        marker2_files: List of marker files from second species
        mapping_file: Gene to protein mapping file
        output_dir: Directory for output files
        gene_prefix: Prefix for gene IDs (default: 'LOC111')
        min_log2fc: Minimum log2 fold change threshold
        max_pval: Maximum p-value threshold
    """
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up logging
    logger = setup_logging(output_dir)
    
    logger.info(f"Starting comparison analysis...")
    logger.info(f"Parameters: gene_prefix={gene_prefix}, min_log2fc={min_log2fc}, max_pval={max_pval}")
    
    # Read mapping file
    logger.info(f"Reading mapping file: {mapping_file}")
    mapping_df = read_file_with_encodings(mapping_file, sep='\t')
    gene_to_protein = dict(zip(mapping_df['gene_id'], mapping_df['protein_id']))
    gene_to_protein = {k: v if v != 'No corresponding protein_id' else None 
                      for k, v in gene_to_protein.items()}
    
    # Store processing results for each file
    marker1_results = {}
    marker2_results = {}
    
    # Process all marker1 files
    for file in marker1_files:
        logger.info(f"Processing marker1 file: {file}")
        proteins = process_marker_file1(file, gene_to_protein, marker1_files,
                                      gene_prefix=gene_prefix,
                                      min_log2fc=min_log2fc,
                                      max_pval=max_pval)
        marker1_results[file] = set(proteins)
        logger.info(f"Found {len(proteins)} specific protein IDs")
    
    # Process all marker2 files
    for file in marker2_files:
        logger.info(f"Processing marker2 file: {file}")
        proteins = process_marker_file2(file, marker2_files)
        marker2_results[file] = set(proteins)
        logger.info(f"Found {len(proteins)} specific protein IDs")
    
    # Compare each pair of files
    comparisons = []
    for marker1_file, marker2_file in product(marker1_files, marker2_files):
        marker1_name = Path(marker1_file).stem
        marker2_name = Path(marker2_file).stem
        
        logger.info(f"\nComparing {marker1_name} and {marker2_name}")
        
        proteins1 = marker1_results[marker1_file]
        proteins2 = marker2_results[marker2_file]
        overlapping = proteins1 & proteins2
        overlap_ratio = len(overlapping) / len(proteins1) if proteins1 else 0
        
        comparison = {
            'marker1_file': marker1_name,
            'marker2_file': marker2_name,
            'marker1_proteins': proteins1,
            'marker2_proteins': proteins2,
            'overlapping_proteins': overlapping,
            'overlap_ratio': overlap_ratio
        }
        
        # Save results for this comparison
        result_file = save_comparison_results(comparison, output_dir)
        logger.info(f"Results saved to: {result_file}")
        
        # Log statistics
        logger.info(f"Number of valid protein IDs in first file: {len(proteins1)}")
        logger.info(f"Number of protein IDs in second file: {len(proteins2)}")
        logger.info(f"Number of overlapping protein IDs: {len(overlapping)}")
        logger.info(f"Overlap ratio: {overlap_ratio:.2%}")
        
        comparisons.append(comparison)
    
    logger.info("Analysis completed!")
    return comparisons

# 使用示例
if __name__ == "__main__":
    # 定义文件路径
    marker1_files = [
        r"D:\nextcloud\pd论文\data\cluster-marker\cluster0_markers.tsv",
        r"D:\nextcloud\pd论文\data\cluster-marker\cluster1_markers.tsv",
        r"D:\nextcloud\pd论文\data\cluster-marker\cluster2_markers.tsv",
        r"D:\nextcloud\pd论文\data\cluster-marker\cluster3_markers.tsv",
        r"D:\nextcloud\pd论文\data\cluster-marker\cluster6_markers.tsv",
        r"D:\nextcloud\pd论文\data\cluster-marker\cluster7_markers.tsv",
    ]
    
    marker2_files = [
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\cnidocyte-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\digestive_filaments-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\neuron-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\mitotic-host-cells-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\alga-hosting-cells-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\calicoblast-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\epidermis-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\gastrodermis-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\germline-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\gland-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\immune-markers.tsv",
        r"D:\nextcloud\pd论文\data\Cell-cluster-marker\unknown-markers.tsv"
    ]
    mapping_file = r"D:\nextcloud\pd论文\data\cluster-marker\loc111_gene_mapping.tsv"
    output_dir = r"D:\nextcloud\pd论文\result\comparison-result2"

    results = compare_marker_pairs(marker1_files, marker2_files, mapping_file, output_dir)