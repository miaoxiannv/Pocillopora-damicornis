#读取gen_list.tsv文件，存成一个列表
import pandas as pd
"""
loc111_genes.tsv文件是LOC111基因的列表
格式如下
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
# 读取loc111_genes.tsv文件
loc111_file = r"D:\nextcloud\pd论文\data\cluster-marker\loc111_genes.tsv"
with open(loc111_file, 'r') as f:
    loc111_genes = [line.strip() for line in f]

print(f"已读取 {len(loc111_genes)} 个LOC111基因")

# 存储gene_id到protein_id的映射
gene_to_protein = {}

"""
gtf文件是基因组注释文件，格式如下
#gtf-version 2.2
#!genome-build Stylophora pistillata v1.1
#!genome-build-accession NCBI_Assembly:GCF_002571385.2
#!annotation-source NCBI Stylophora pistillata Annotation Release 100
NW_019217784.1	Gnomon	gene	30	5662	.	-	.	gene_id "LOC111326392"; transcript_id ""; db_xref "GeneID:111326392"; gbkey "Gene"; gene "LOC111326392"; gene_biotype "lncRNA"; 
NW_019217784.1	Gnomon	transcript	30	5662	.	-	.	gene_id "LOC111326392"; transcript_id "XR_002694437.1"; db_xref "GeneID:111326392"; gbkey "ncRNA"; gene "LOC111326392"; model_evidence "Supporting evidence includes similarity to: 100% coverage of the annotated genomic feature by RNAseq alignments, including 6 samples with support for all annotated introns"; product "uncharacterized LOC111326392"; transcript_biotype "lnc_RNA"; 
NW_019217784.1	Gnomon	exon	5214	5662	.	-	.	gene_id "LOC111326392"; transcript_id "XR_002694437.1"; db_xref "GeneID:111326392"; gene "LOC111326392"; model_evidence "Supporting evidence includes similarity to: 100% coverage of the annotated genomic feature by RNAseq alignments, including 6 samples with support for all annotated introns"; product "uncharacterized LOC111326392"; transcript_biotype "lnc_RNA"; exon_number "1"; 
NW_019217784.1	Gnomon	exon	4427	4584	.	-	.	gene_id "LOC111326392"; transcript_id "XR_002694437.1"; db_xref "GeneID:111326392"; gene "LOC111326392"; model_evidence "Supporting evidence includes similarity to: 100% coverage of the annotated genomic feature by RNAseq alignments, including 6 samples with support for all annotated introns"; product "uncharacterized LOC111326392"; transcript_biotype "lnc_RNA"; exon_number "2"; 
NW_019217784.1	Gnomon	exon	30	1561	.	-	.	gene_id "LOC111326392"; transcript_id "XR_002694437.1"; db_xref "GeneID:111326392"; gene "LOC111326392"; model_evidence "Supporting evidence includes similarity to: 100% coverage of the annotated genomic feature by RNAseq alignments, including 6 samples with support for all annotated introns"; product "uncharacterized LOC111326392"; transcript_biotype "lnc_RNA"; exon_number "3"; 
NW_019217784.1	Gnomon	gene	2759	16078	.	+	.	gene_id "LOC111326177"; transcript_id ""; db_xref "GeneID:111326177"; gbkey "Gene"; gene "LOC111326177"; gene_biotype "protein_coding"; 
NW_019217784.1	Gnomon	transcript	2759	16078	.	+	.	gene_id "LOC111326177"; transcript_id "XM_022930150.1"; db_xref "GeneID:111326177"; gbkey "mRNA"; gene "LOC111326177"; model_evidence "Supporting evidence includes similarity to: 1 EST, 2 Proteins, and 100% coverage of the annotated genomic feature by RNAseq alignments, including 8 samples with support for all annotated introns"; product "heat shock factor protein 1-like, transcript variant X1"; transcript_biotype "mRNA"; 
NW_019217784.1	Gnomon	exon	2759	2893	.	+	.	gene_id "LOC111326177"; transcript_id "XM_022930150.1"; db_xref "GeneID:111326177"; gene "LOC111326177"; model_evidence "Supporting evidence includes similarity to: 1 EST, 2 Proteins, and 100% coverage of the annotated genomic feature by RNAseq alignments, including 8 samples with support for all annotated introns"; product "heat shock factor protein 1-like, transcript variant X1"; transcript_biotype "mRNA"; exon_number "1"; 
NW_019217784.1	Gnomon	exon	7187	7359	.	+	.	gene_id "LOC111326177"; transcript_id "XM_022930150.1"; db_xref "GeneID:111326177"; gene "LOC111326177"; model_evidence "Supporting evidence includes similarity to: 1 EST, 2 Proteins, and 100% coverage of the annotated genomic feature by RNAseq alignments, including 8 samples with support for all annotated introns"; product "heat shock factor protein 1-like, transcript variant X1"; transcript_biotype "mRNA"; exon_number "2"; 
NW_019217784.1	Gnomon	exon	7521	7629	.	+	.	gene_id "LOC111326177"; transcript_id "XM_022930150.1"; db_xref "GeneID:111326177"; gene "LOC111326177"; model_evidence "Supporting evidence includes similarity to: 1 EST, 2 Proteins, and 100% coverage of the annotated genomic feature by RNAseq alignments, including 8 samples with support for all annotated introns"; product "heat shock factor protein 1-like, transcript variant X1"; transcript_biotype "mRNA"; exon_number "3"; 
NW_019217784.1	Gnomon	CDS	10573	10631	.	+	0	gene_id "LOC111326177"; transcript_id "XM_022930150.1"; db_xref "GeneID:111326177"; gbkey "CDS"; gene "LOC111326177"; product "heat shock factor protein-like isoform X1"; protein_id "XP_022785885.1"; exon_number "7"; 
NW_019217784.1	Gnomon	CDS	10978	11092	.	+	1	gene_id "LOC111326177"; transcript_id "XM_022930150.1"; db_xref "GeneID:111326177"; gbkey "CDS"; gene "LOC111326177"; product "heat shock factor protein-like isoform X1"; protein_id "XP_022785885.1"; exon_number "8"; 
NW_019217784.1	Gnomon	CDS	11708	11768	.	+	0	gene_id "LOC111326177"; transcript_id "XM_022930150.1"; db_xref "GeneID:111326177"; gbkey "CDS"; gene "LOC111326177"; product "heat shock factor protein-like isoform X1"; protein_id "XP_022785885.1"; exon_number "9"; 

"""

# 读取gtf文件并查找对应的protein_id
gtf_file = r'D:\nextcloud\pd论文\data\NCBI-SP\GCF_002571385.2_Stylophora_pistillata_v1.1_genomic.gtf\GCF_0025.2_S'
with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):  # 跳过注释行
            continue
        if '\tCDS\t' not in line:  # 只处理CDS行
            continue
            
        # 解析gene_id和protein_id
        attributes = line.strip().split('\t')[-1]
        gene_id = None
        protein_id = None
        
        # 解析属性字段
        for attr in attributes.split('; '):
            if 'gene_id' in attr:
                gene_id = attr.split('"')[1]
            elif 'protein_id' in attr:
                protein_id = attr.split('"')[1]
                
        if gene_id and protein_id:
            gene_to_protein[gene_id] = protein_id

# 输出结果到tsv文件
output_file = r'D:\nextcloud\pd论文\data\cluster-marker\loc111_gene_mapping.tsv'
with open(output_file, 'w') as out:
    # 写入表头
    out.write("gene_id\tprotein_id\n")
    # 写入数据
    for gene in loc111_genes:
        if gene in gene_to_protein:
            out.write(f"{gene}\t{gene_to_protein[gene]}\n")
        else:
            out.write(f"{gene}\t未找到对应的protein_id\n")

# 打印统计信息
found_count = sum(1 for gene in loc111_genes if gene in gene_to_protein)
print(f"总共有 {len(loc111_genes)} 个LOC111基因")
print(f"找到对应protein_id的基因数量：{found_count}")
print(f"未找到对应protein_id的基因数量：{len(loc111_genes) - found_count}")
print(f"结果已保存到文件：{output_file}")
