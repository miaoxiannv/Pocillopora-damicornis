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


def read_gene_list(file_path):
    """Read gene list from file"""
    with open(file_path, 'r') as f:
        return [line.strip() for line in f]


def parse_gtf_for_protein_ids(gtf_file):
    """Parse GTF file and extract gene to protein ID mappings"""
    gene_to_protein = {}
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if '\tCDS\t' not in line:
                continue
            
            attributes = line.strip().split('\t')[-1]
            gene_id = None
            protein_id = None
            
            for attr in attributes.split('; '):
                if 'gene_id' in attr:
                    gene_id = attr.split('"')[1]
                elif 'protein_id' in attr:
                    protein_id = attr.split('"')[1]
            
            if gene_id and protein_id:
                gene_to_protein[gene_id] = protein_id
    
    return gene_to_protein


def write_mapping_results(output_file, loc111_genes, gene_to_protein):
    """Write gene to protein ID mapping results to file"""
    with open(output_file, 'w') as out:
        out.write("gene_id\tprotein_id\n")
        for gene in loc111_genes:
            if gene in gene_to_protein:
                out.write(f"{gene}\t{gene_to_protein[gene]}\n")
            else:
                out.write(f"{gene}\tNo corresponding protein_id found\n")


if __name__ == "__main__":
    loc111_file = r"D:\nextcloud\pd论文\data\merge-data\SP-NCBI-data\genlist.tsv"
    gtf_file = r'D:\nextcloud\pd论文\data\NCBI-SP\GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.gtf\GCF_0025.2_S'
    output_file = r'D:\nextcloud\pd论文\data\cluster-marker\ncbi_gene_mapping.tsv'
    
    loc111_genes = read_gene_list(loc111_file)
    print(f"Read {len(loc111_genes)} LOC111 genes")
    
    gene_to_protein = parse_gtf_for_protein_ids(gtf_file)
    
    write_mapping_results(output_file, loc111_genes, gene_to_protein)
    
    found_count = sum(1 for gene in loc111_genes if gene in gene_to_protein)
    print(f"Total LOC111 genes: {len(loc111_genes)}")
    print(f"Number of genes with protein_id: {found_count}")
    print(f"Number of genes without protein_id: {len(loc111_genes) - found_count}")
    print(f"Results saved to file: {output_file}")
