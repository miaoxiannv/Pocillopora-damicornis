"""
Gene to Protein ID Mapper.
This script maps gene IDs (LOC111 format) to their corresponding protein IDs using GTF annotation file.

Input files:
    - gene_list.tsv: List of LOC111 format gene IDs
    - GTF annotation file: Contains gene and protein information

Output:
    - ncbi_gene_mapping.tsv: Mapping between gene IDs and protein IDs

Usage:
    python gene_to_protein.py

Author: Shengyao Zhang
Date: 2024-12-19
Version: 1.0
"""

import pandas as pd

"""
Example of loc111_genes.tsv format:
LOC111325671
LOC111326254
LOC111326848
LOC111339103
LOC111343026
LOC111322693
LOC111320871
LOC111340891
LOC111327082
LOC111320821
"""
# Read loc111_genes.tsv file
loc111_file = r"D:\nextcloud\pd论文\data\merge-data\SP-NCBI-data\genlist.tsv"
with open(loc111_file, 'r') as f:
    loc111_genes = [line.strip() for line in f]

print(f"Read {len(loc111_genes)} LOC111 genes")

# Store gene_id to protein_id mapping
gene_to_protein = {}

"""
GTF file format example:
#gtf-version 2.2
#!genome-build Stylophora pistillata v1.1
#!genome-build-accession NCBI_Assembly:GCF_002571385.2
#!annotation-source NCBI Stylophora pistillata Annotation Release 100
NW_019217784.1	Gnomon	gene	30	5662	.	-	.	gene_id "LOC111326392"; ...
"""

# Read GTF file and find corresponding protein_ids
gtf_file = r'D:\nextcloud\pd论文\data\NCBI-SP\GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.gtf\GCF_0025.2_S'
with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):  # Skip comment lines
            continue
        if '\tCDS\t' not in line:  # Only process CDS lines
            continue
            
        # Parse gene_id and protein_id
        attributes = line.strip().split('\t')[-1]
        gene_id = None
        protein_id = None
        
        # Parse attribute fields
        for attr in attributes.split('; '):
            if 'gene_id' in attr:
                gene_id = attr.split('"')[1]
            elif 'protein_id' in attr:
                protein_id = attr.split('"')[1]
                
        if gene_id and protein_id:
            gene_to_protein[gene_id] = protein_id

# Output results to TSV file
output_file = r'D:\nextcloud\pd论文\data\cluster-marker\ncbi_gene_mapping.tsv'
with open(output_file, 'w') as out:
    # Write header
    out.write("gene_id\tprotein_id\n")
    # Write data
    for gene in loc111_genes:
        if gene in gene_to_protein:
            out.write(f"{gene}\t{gene_to_protein[gene]}\n")
        else:
            out.write(f"{gene}\tNo corresponding protein_id found\n")

# Print statistics
found_count = sum(1 for gene in loc111_genes if gene in gene_to_protein)
print(f"Total LOC111 genes: {len(loc111_genes)}")
print(f"Number of genes with protein_id: {found_count}")
print(f"Number of genes without protein_id: {len(loc111_genes) - found_count}")
print(f"Results saved to file: {output_file}")
