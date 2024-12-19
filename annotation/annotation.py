"""
Annotation module for processing gene annotations from various databases.
This script provides base classes and implementations for handling different types of annotations
including GO, KEGG, Pfam, and KOG annotations.

Usage:
    Import this module to use annotation functionality:
    from annotation import AnnotationGO, AnnotationKegg, AnnotationPfam, AnnotationKog

Author: Shengyao Zhang
Date: 2024-12-19
Version: 1.0
"""

import os
import pandas as pd
# Import abstract method
from abc import ABC, abstractmethod

class GeneList(object):
    def __init__(self, gene_list_file):
        self.gene_list_file = gene_list_file
        self.gene_list = None
        self.read_gene_list()

    def read_gene_list(self):
        gene_df = pd.read_csv(self.gene_list_file, header=None, index_col=False)
        gene_df.columns = ['gene']
        # Lower the gene column
        gene_df['gene'] = gene_df['gene'].str.strip().str.lower()
        # Transform the gene column to python list
        self.gene_list = gene_df['gene'].tolist()
    
    def get_gene_list(self):
        return self.gene_list

class Annotation(ABC):
    def __init__(self, anno_file):
        self.anno_file = anno_file
        self.anno_df = None
    
    # Read annotation file
    def read_anno_file(self):
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
        pass
    
    def search(self, gene_list):
        # Convert gene names to lowercase for matching
        gene_list = [gene.strip().lower() for gene in gene_list]
        
        # Use reindex to find annotation information
        # reindex will automatically handle non-existent genes by setting them to NaN
        res = self.anno_df.reindex(gene_list)
        found_genes = res.dropna(how='all').shape[0]
        # Reset index
        res.reset_index(inplace=True)
        # Count found and not found genes
        
        not_found_genes = res.shape[0] - found_genes
        
        print(f"Search Results:")
        print(f"- Total input genes: {len(gene_list)}")
        print(f"- Genes with annotation: {found_genes}")
        print(f"- Genes without annotation: {not_found_genes}")
        
        return res
    
    # Export the result to excel or tsv file
    def export(self, res, out_file):
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

class AnnotationGO(Annotation):
    '''
    the annotation file format:
    Gene_id	GO_number	GO_id; GO_description; GO_class
    transcript_HQ_P2_transcript0/f8p0/9897	1	GO:0005515; protein binding; molecular_function
    transcript_HQ_P2_transcript1/f2p0/9329	10	GO:0006914; autophagy; biological_process
    transcript_HQ_P2_transcript10/f2p0/8399	4	GO:0004222; metalloendopeptidase activity; molecular_function
    transcript_HQ_P2_transcript100/f5p0/6808	8	GO:0016192; vesicle-mediated transport; biological_process
    transcript_HQ_P2_transcript1000/f3p0/4849	7	GO:0016787; hydrolase activity; molecular_function
    transcript_HQ_P2_transcript10000/f12p0/2537	4	GO:0005666; DNA-directed RNA polymerase III complex; cellular_component
    transcript_HQ_P2_transcript10001/f2p0/2514	4	GO:0015031; protein transport; biological_process
    '''
    def __init__(self, anno_file):
        super().__init__(anno_file)

    def set_index(self):
        #strip the gene_id and lower
        self.anno_df['Gene_id'] = self.anno_df['Gene_id'].str.strip().str.lower()
        self.anno_df.set_index('Gene_id', inplace=True)
        return self.anno_df

class AnnotationKegg(Annotation):
    '''
    the annotation file format:
    Query_id	Subject_id	KO_ID	KO_NAME	KO_DEFINITION	KO_EC	KO_PATHWAY
    transcript_HQ_P2_transcript10902/f3p0/2417	nve:NEMVE_v1g112567	K00830	AGXT	alanine-glyoxylate transaminase / serine-glyoxylate transaminase / serine-pyruvate transaminase [EC:2.6.1.44 2.6.1.45 2.6.1.51]	2.6.1.44,2.6.1.45,2.6.1.51	ko01200; Metabolism; Global and overview maps; Carbon metabolism | ko00630; Metabolism; Carbohydrate metabolism; Glyoxylate and dicarboxylate metabolism | ko00680; Metabolism; Energy metabolism; Methane metabolism | ko00250; Metabolism; Amino acid metabolism; Alanine, aspartate and glutamate metabolism | ko00260; Metabolism; Amino acid metabolism; Glycine, serine and threonine metabolism | ko04146; Cellular Processes; Transport and catabolism; Peroxisome
    transcript_HQ_P2_transcript12712/f6p0/2221	cfr:102522288	-	-	-	-	-
    transcript_HQ_P2_transcript14522/f12p0/2021	nve:NEMVE_v1g173936	K00552	GNMT	glycine N-methyltransferase	2.1.1.20	ko00260; Metabolism; Amino acid metabolism; Glycine, serine and threonine metabolism
    '''

    def __init__(self, anno_file):
        super().__init__(anno_file) 

    def set_index(self):
        self.anno_df['Query_id'] = self.anno_df['Query_id'].str.strip().str.lower()
        
        # Check for duplicates
        duplicates = self.anno_df['Query_id'].duplicated().sum()
        print(f"Found {duplicates} duplicate gene IDs")
        
        # Group and merge duplicates
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
        return self.anno_df

class AnnotationPfam(Annotation):
    '''
    Gene_id	Pfam_number	Pfam_id:Pfam_description
    transcript_HQ_P2_transcript5854/f2p0/3088	3	PF01344:Kelch motif
    transcript_HQ_P2_transcript2813/f7p0/3807	1	PF00400:WD domain, G-beta repeat
    transcript_HQ_P2_transcript18931/f16p0/1250	2	PF02045:CCAAT-binding transcription factor (CBF-B/NF-YA) subunit B
    transcript_HQ_P2_transcript2528/f2p0/3909	1	PF00641:Zn-finger in Ran binding protein and others
    transcript_HQ_P2_transcript5639/f3p0/3108	2	PF00069:Protein kinase domain
    transcript_HQ_P2_transcript2186/f2p0/4083	4	PF00071:Ras family
    transcript_HQ_P2_transcript6823/f2p0/2906	1	PF00620:RhoGAP domain
    transcript_HQ_P2_transcript962/f2p0/4890	4	PF00682:HMGL-like
    transcript_HQ_P2_transcript18648/f2p0/1352	5	PF00292:'Paired box' domain
    transcript_HQ_P2_transcript17810/f2p0/1552	4	PF00095:WAP-type (Whey Acidic Protein) 'four-disulfide core'
    transcript_HQ_P2_transcript13051/f2p0/2188	1	PF05887:Procyclic acidic repetitive protein (PARP)
    transcript_HQ_P2_transcript12396/f6p0/2239	3	PF04088:Peroxin 13, N-terminal region
    '''
    def __init__(self, anno_file):
        super().__init__(anno_file)

    def set_index(self):
        #strip the gene_id and lower
        self.anno_df['Gene_id'] = self.anno_df['Gene_id'].str.strip().str.lower()
        self.anno_df.set_index('Gene_id', inplace=True)
        return self.anno_df

class AnnotationKog(Annotation):
    '''
    Gene_id	Identity	E_value	KOG_gene_id	KOG_num	Functional_description	Functional_class	Class_description
    transcript_HQ_P2_transcript10902/f3p0/2417	48.7	6.10E-113	Hs4557289	KOG2862	Alanine-glyoxylate aminotransferase AGT1	R	General function prediction only ;
    transcript_HQ_P2_transcript12713/f5p0/2200	49.6	2.30E-21	Hs22041557	KOG2723	Uncharacterized conserved protein, contains BTB/POZ domain	R	General function prediction only ;
    transcript_HQ_P2_transcript14523/f2p0/2027	41.2	4.80E-26	7301068	KOG3335	Predicted coiled-coil protein	R	General function prediction only ;
    transcript_HQ_P2_transcript13619/f3p0/2127	80	2.20E-66	Hs20560164	KOG2157	Predicted tubulin-tyrosine ligase	O	Posttranslational modification, protein turnover, chaperones ;
    transcript_HQ_P2_transcript10905/f2p0/2440	45.8	1.10E-125	Hs19263340	KOG0258	Alanine aminotransferase	E	Amino acid transport and metabolism ;

    '''
    def __init__(self, anno_file):
        super().__init__(anno_file)
    
    def set_index(self):
        #strip the gene_id and lower
        self.anno_df['Gene_id'] = self.anno_df['Gene_id'].str.strip().str.lower()
        self.anno_df.set_index('Gene_id', inplace=True)
        return self.anno_df


if __name__ == "__main__":
    # Test file path
    gene_list_file = r"C:\Users\伯言\Desktop\test001\filtered_genes.txt"
    go_anno_file = r"D:\nextcloud\pd论文\data\test\annotation\GO.anno.xls"
    kegg_anno_file = r"D:\nextcloud\pd论文\data\test\annotation\P2.KEGG.filter.m8.anno.xls"
    pfam_anno_file = r"D:\nextcloud\pd论文\data\test\annotation\P2.pfam.anno.xls"
    kog_anno_file = r"D:\nextcloud\pd论文\data\test\annotation\P2.KOG.filter.m8.anno.xls"
    
    # Create output directory
    out_dir = r"D:\nextcloud\pd论文\result\test\annotation"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Create gene list object
    gene_list = GeneList(gene_list_file)
    gene_list = gene_list.get_gene_list()

    # Test GO annotation
    print("Testing GO annotation...")
    go_anno = AnnotationGO(go_anno_file)
    go_anno.read_anno_file()
    go_anno.set_index()
    go_results = go_anno.search(gene_list)
    print(f"GO annotation result row count: {len(go_results)}")
    print(go_results.head())
    go_anno.export(go_results, os.path.join(out_dir, "go_results.tsv"))
    print("\n" + "="*50 + "\n")

    # Test KEGG annotation
    print("Testing KEGG annotation...")
    kegg_anno = AnnotationKegg(kegg_anno_file)
    kegg_anno.read_anno_file()
    kegg_anno.set_index()
    kegg_results = kegg_anno.search(gene_list)
    print(f"KEGG annotation result row count: {len(kegg_results)}")
    print(kegg_results.head())
    kegg_anno.export(kegg_results, os.path.join(out_dir, "kegg_results.tsv"))
    print("\n" + "="*50 + "\n")

    # Test Pfam annotation
    print("Testing Pfam annotation...")
    pfam_anno = AnnotationPfam(pfam_anno_file)
    pfam_anno.read_anno_file()
    pfam_anno.set_index()
    pfam_results = pfam_anno.search(gene_list)
    print(f"Pfam annotation result row count: {len(pfam_results)}")
    print(pfam_results.head())
    pfam_anno.export(pfam_results, os.path.join(out_dir, "pfam_results.tsv"))
    print("\n" + "="*50 + "\n")

    # Test KOG annotation
    print("Testing KOG annotation...")
    kog_anno = AnnotationKog(kog_anno_file)
    kog_anno.read_anno_file()
    kog_anno.set_index()
    kog_results = kog_anno.search(gene_list)
    print(f"KOG annotation result row count: {len(kog_results)}")
    print(kog_results.head())
    kog_anno.export(kog_results, os.path.join(out_dir, "kog_results.tsv"))
    print("\n" + "="*50 + "\n")
























