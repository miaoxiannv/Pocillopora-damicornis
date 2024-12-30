# Bioinformatics Analysis Pipeline

This project is a bioinformatics analysis pipeline that includes multiple steps from generating reference genome, processing 10X genomics data, gene annotation to cell type identification. Here are detailed descriptions for each part.

## 1. Generating GTF Files

The `full_length_make_gtf_ref.py` script converts FASTA files to GTF and new FASTA format. It takes a FASTA file or a directory containing FASTA files as input, and generates a new FASTA file and a GTF file based on the input sequences.

### Features

- Supports single FASTA file or directory containing multiple FASTA files as input
- Generates a new FASTA file with cleaned up sequence names
- Creates a GTF file with gene and transcript annotations
- Handles duplicate sequence names by appending a unique index
- Provides command-line arguments for specifying input and output paths

### Usage
python full_length_make_gtf_ref.py -i <input_fasta> -o <output_directory>

## 2. 10X Genomics Data Processing

This Python script merges 10X Genomics data from different species into a unified format. It processes barcode, feature, and matrix files from multiple species and combines them into a single dataset.

### Features

- Supports merging data from multiple species
- Handles input files in plain text or gzip compressed format
- Generates unified output files: barcodes.tsv, features.tsv, and matrix.mtx
- Provides error handling and logging for data processing issues

### Usage

1. Prepare the input data files for each species:
   - barcode.tsv: Contains cell barcodes
   - feature.tsv: Contains gene/feature information
   - matrix.mtx: Contains expression matrix in sparse format
2. Update the data_paths dictionary in the script with the file paths for each species

## 3. Gene Annotation 

This Python script annotates marker genes with additional information from various databases such as GO, KEGG, Pfam, and KOG. It reads a marker gene file and multiple annotation files, merges the information, and exports the annotated marker genes to output files.

### Features

- Supports annotation from multiple databases: GO, KEGG, Pfam, and KOG
- Reads marker gene file and annotation files in TSV format
- Merges annotation information with marker gene data
- Handles duplicate gene IDs and merges their annotations
- Exports annotated results to TSV files
- Provides summary statistics of annotation results

### Usage

1. Prepare the input files:
   - Marker gene file in TSV format
   - Annotation files in TSV format for GO, KEGG, Pfam, and KOG
2. Update the file paths in the `__main__` section of the script

## 4. Cell Type Identification and Gene Name Processing

This part includes three Python scripts for cell type identification and gene name processing:

1. `celltype_identification2.py`: Compares marker gene files from different species to identify overlapping proteins.
2. `gene_to_protein.py`: Maps gene IDs (LOC111 format) to corresponding protein IDs using a GTF annotation file.  
3. `extract_gene_names.py`: Extracts gene names from marker gene files and saves LOC111 format gene names separately.

### Features

- Compares marker gene files from different species to identify overlapping proteins
- Maps LOC111 format gene IDs to corresponding protein IDs
- Extracts gene names from marker gene files and saves LOC111 format gene names separately

### Usage 

1. Prepare input files:
   - Marker gene files from different species (tsv format)
   - Gene ID to protein ID mapping file (tsv format)
   - GTF annotation file
2. Update file paths in the scripts
3. Run the scripts 
4. Check the output results

## 5. GSEA Enrichment Analysis

This is a Python module for performing Gene Set Enrichment Analysis (GSEA).

### Features  

- Perform hypergeometric test on gene lists and gene sets
- P-value correction for multiple hypothesis testing
- Support multiple gene set database formats such as GMT, GO, KEGG, KOG, and Pfam
- Flexible input and output options
Here's the code with the markdown format corrected:

### Usage Example
```python
from gsea.algorithms.hypergeom import calcu_hypergeom
from gsea.data.gene_list_obj import GeneList_Obj
from gsea.data.gene_set_obj_kegg import Gmt_stat

# Load gene list and gene set
gene_list = GeneList_Obj("path/to/gene_list.txt").gene_list
gmt_obj = Gmt_stat.from_gmt("path/to/gene_set.gmt")

# Perform GSEA analysis
res = calcu_hypergeom(gene_list, gmt_obj, min_count=5)

# Output enrichment results
print(res)
```

### Dependencies

- Python 3.6+
- NumPy
- Pandas
- SciPy
- statsmodels
# 生物信息学分析流程

本项目是一个生物信息学分析流程,包含从生成参考基因组、处理10X基因组学数据、基因注释到细胞类型鉴定等多个步骤。以下是每个部分的详细说明。

## 1. 生成GTF文件

`full_length_make_gtf_ref.py`脚本将FASTA文件转换为GTF和新的FASTA格式。它支持单个FASTA文件或包含多个FASTA文件的目录作为输入,并基于输入序列生成新的FASTA文件和GTF文件。

### 功能

- 支持单个FASTA文件或包含多个FASTA文件的目录作为输入 
- 生成带有清理后序列名称的新FASTA文件
- 创建包含基因和转录本注释的GTF文件
- 通过附加唯一索引来处理重复的序列名称
- 提供命令行参数以指定输入和输出路径

### 使用方法
python full_length_make_gtf_ref.py -i <input_fasta> -o <output_directory>

## 2. 10X Genomics数据处理

此Python脚本将来自不同物种的10X Genomics数据合并为统一格式。它处理来自多个物种的barcode、feature和matrix文件,并将它们组合成一个统一的数据集。 

### 功能

- 支持合并来自多个物种的数据
- 处理纯文本或gzip压缩格式的输入文件
- 生成统一的输出文件:barcodes.tsv、features.tsv和matrix.mtx
- 提供数据处理问题的错误处理和日志记录

### 使用方法

1. 为每个物种准备输入数据文件:
   - barcode.tsv:包含细胞条形码
   - feature.tsv:包含基因/特征信息
   - matrix.mtx:包含稀疏格式的表达矩阵
2. 更新脚本中的data_paths字典,使用每个物种的文件路径

## 3. 基因注释

这个Python脚本使用来自GO、KEGG、Pfam和KOG等多个数据库的额外信息对标记基因进行注释。它读取标记基因文件和多个注释文件,合并信息,并将注释后的标记基因导出到输出文件。

### 功能

- 支持来自多个数据库的注释:GO、KEGG、Pfam和KOG
- 读取TSV格式的标记基因文件和注释文件
- 将注释信息与标记基因数据合并
- 处理重复的基因ID并合并其注释
- 将注释结果导出到TSV文件
- 提供注释结果的汇总统计信息

### 使用方法

1. 准备输入文件:
   - TSV格式的marker基因文件
   - TSV格式的GO、KEGG、Pfam和KOG注释文件
2. 更新脚本__main__部分中的文件路径

## 4. 细胞类型鉴定和基因名称处理

这部分包括三个Python脚本,用于细胞类型鉴定和基因名称处理:

1. `celltype_identification2.py`:比较不同物种的标记基因文件,识别重叠的蛋白质。
2. `gene_to_protein.py`:使用GTF注释文件将基因ID(LOC111格式)映射到对应的蛋白质ID。

3. `extract_gene_names.py`:从标记基因文件中提取基因名称,并将LOC111格式的基因名称单独保存。


## 5. GSEA富集分析

这是一个用于执行基因集富集分析(GSEA)的Python模块。

### 功能

- 对基因列表和基因集执行超几何检验
- 多重假设检验的p值校正
- 支持多种基因集数据库格式,如GMT、GO、KEGG、KOG和Pfam
- 灵活的输入输出选项

### 使用示例
```python
from gsea.algorithms.hypergeom import calcu_hypergeom
from gsea.data.gene_list_obj import GeneList_Obj
from gsea.data.gene_set_obj_kegg import Gmt_stat

# 加载基因列表和基因集
gene_list = GeneList_Obj("path/to/gene_list.txt").gene_list
gmt_obj = Gmt_stat.from_gmt("path/to/gene_set.gmt")

# 执行GSEA分析
res = calcu_hypergeom(gene_list, gmt_obj, min_count=5)

# 输出富集结果
print(res)
```

### 依赖

- Python 3.6+
- NumPy
- Pandas
- SciPy
- statsmodels