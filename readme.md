# Bioinformatics Analysis Pipeline

This project is a comprehensive bioinformatics analysis pipeline that includes multiple steps: generating reference genomes, processing 10X Genomics data, gene annotation, and cell type identification. Below are detailed descriptions for each component.

## 1. Generating GTF Files

The `full_length_make_gtf_ref.py` script converts FASTA files to GTF and new FASTA format. It takes a FASTA file or a directory containing FASTA files as input and generates a new FASTA file and a GTF file based on the input sequences.

### Features

- Supports single FASTA file or directory containing multiple FASTA files as input
- Generates a new FASTA file with cleaned up sequence names
- Creates a GTF file with gene and transcript annotations
- Handles duplicate sequence names by appending a unique index
- Provides command-line arguments for specifying input and output paths

### Usage

```bash
python full_length_make_gtf_ref.py -i <input_fasta> -o <output_directory>
```

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
- Reads marker gene files and annotation files in TSV format
- Merges annotation information with marker gene data
- Handles duplicate gene IDs and merges their annotations
- Exports annotated results to TSV files
- Provides summary statistics of annotation results

### Usage

1. Prepare the input files:
   - Marker gene file in TSV format
   - Annotation files in TSV format for GO, KEGG, Pfam, and KOG
2. Update the file paths in the `__main__` section of the script
3. Run the script:

```bash
python annotation_with_marker_info.py
```

## 4. Cell Type Identification and Gene Name Processing

This part includes three Python scripts for cell type identification and gene name processing:

1. `celltype_identification2.py`: Compares marker gene files from different species to identify overlapping proteins
2. `gene_to_protein.py`: Maps gene IDs (LOC111 format) to corresponding protein IDs using a GTF annotation file
3. `extract_gene_names.py`: Extracts gene names from marker gene files and saves LOC111 format gene names separately

### Features

- Compares marker gene files from different species to identify overlapping proteins
- Maps LOC111 format gene IDs to corresponding protein IDs
- Extracts gene names from marker gene files and saves LOC111 format gene names separately

### Usage

1. Prepare input files:
   - Marker gene files from different species (TSV format)
   - Gene ID to protein ID mapping file (TSV format)
   - GTF annotation file
2. Update file paths in the scripts
3. Run the scripts:

```bash
python celltype_identification2.py
python gene_to_protein.py
python extract_gene_names.py
```

4. Check the output results

## 5. GSEA Enrichment Analysis

This is a Python module for performing Gene Set Enrichment Analysis (GSEA).

### Features

- Perform hypergeometric test on gene lists and gene sets
- P-value correction for multiple hypothesis testing
- Support multiple gene set database formats such as GMT, GO, KEGG, KOG, and Pfam
- Flexible input and output options

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

## Installation

1. Clone the repository:

```bash
git clone <repository-url>
cd pd
```

2. Install dependencies:

```bash
pip install -r requirement.txt
```

## Project Structure

```
.
├── align_count/              # GTF and FASTA file generation
├── annotation/               # Gene annotation tools
├── data_processing/          # 10X Genomics data processing
├── gsea/                     # GSEA analysis module
│   ├── algorithms/          # Analysis algorithms
│   └── data/                # Data objects and loaders
├── processing_gene_names/    # Gene name processing utilities
├── readme.md                 # This file
└── requirement.txt          # Python dependencies
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

[Add your license information here]

## Contact

[Add your contact information here]
