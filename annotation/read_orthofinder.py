"""
OrthoFinder results processor.
This script reads OrthoFinder output files, extracts gene names from orthogroups,
and exports them to a simplified format.

Usage:
    python read_orthofinder.py
    This will process the specified OrthoFinder output file and export unique gene names.

Author: Shengyao Zhang
Date: 2024-12-19
Version: 1.0
"""

def read_orthogroup_genes(file_path):
    """Read orthogroup file and extract gene names from the third column into a dictionary"""
    
    # Read file using tab delimiter, skip header
    with open(file_path, 'r') as f:
        next(f)
        lines = f.readlines()
    
    gene_dict = {}
    for line in lines:
        # Skip empty lines
        if not line.strip():
            continue
            
        # Split line by tab
        columns = line.strip().split('\t')
        
        # Ensure at least 3 columns exist
        if len(columns) >= 3:
            # Get genes from third column
            genes_str = columns[2]
            if genes_str:  # Ensure not empty string
                # Split gene names by comma and space
                genes = genes_str.split(',')
                # Clean each gene name and store in dictionary
                for gene in genes:
                    gene = gene.strip()
                    if gene:  # Ensure gene name is not empty
                        gene_dict[gene] = True
    
    print(f"Found {len(gene_dict)} unique genes in total")
    
    return gene_dict

def export_genes(gene_dict, output_file):
    """Export gene dictionary to file, one gene per line"""
    with open(output_file, 'w') as f:
        for gene in sorted(gene_dict.keys()):
            f.write(f"{gene}\n")
    print(f"Gene list has been saved to: {output_file}")

if __name__ == "__main__":
    input_file = r"D:\nextcloud\pd论文\data\test\Cell_type_prediction\PD-SP-cell\PD-SP-cell_Orthogroups.tsv"  # Replace with your input file path
    output_file = r"D:\nextcloud\pd论文\data\test\Cell_type_prediction\PD-SP-cell\PD-SP-cell_Orthogroups_genes.txt"  # Replace with your desired output file path
    
    # Read genes
    gene_dict = read_orthogroup_genes(input_file)
    
    # Export results
    if gene_dict:
        export_genes(gene_dict, output_file)