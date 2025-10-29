"""
Extract gene names from marker files.

This script reads marker gene files and extracts unique gene names,
with special handling for LOC111 format gene IDs.

Marker file format example:
"p_val" "avg_log2FC" "pct.1" "pct.2" "p_val_adj" "cluster" "gene"
"LOC131768943" 5.98636196686601e-135 -4.71546249021902 0.016 0.352 3.21234169503997e-130 "0" "LOC131768943"
"LOC113677177" 3.38178484945174e-130 -3.41989240052076 0.02 0.354 1.8146995680643e-125 "0" "LOC113677177"
"transcript-HQ-P2-transcript13875/f192p0/1996" 2.35948385593556e-124 -2.26149823236855 0.025 0.355 1.26612263193358e-119 "0" "transcript-HQ-P2-transcript13875/f192p0/1996"

Author: Shengyao Zhang
Date: 2024-12-19
Version: 1.0
"""
import pandas as pd


if __name__ == "__main__":
    marker_file = r"D:\nextcloud\pd论文\data\cluster-marker\allmarkers.tsv"
    marker_df = pd.read_csv(marker_file, sep=" ")
    print(marker_df.head())
    
    gene_list = []
    loc111_genes = []
    
    for gene in marker_df['gene'].unique():
        gene_clean = gene.strip('"')
        gene_list.append(gene_clean)
        if gene_clean.startswith("LOC111"):
            loc111_genes.append(gene_clean)
    
    output_file = r"D:\nextcloud\pd论文\data\cluster-marker\gene_list.tsv"
    with open(output_file, 'w') as f:
        for gene in gene_list:
            f.write(f"{gene}\n")
    print(f"Extracted {len(gene_list)} unique gene names and saved to {output_file}")
    
    output_file = r"D:\nextcloud\pd论文\data\cluster-marker\loc111_genes.tsv"
    with open(output_file, 'w') as f:
        for gene in loc111_genes:
            f.write(f"{gene}\n")
    print(f"Extracted {len(loc111_genes)} LOC111 gene names and saved to {output_file}")