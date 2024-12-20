"""
Convert FASTA to GTF and new FASTA.
This script takes a FASTA file or a directory containing FASTA files as input,
and generates a new FASTA file and a GTF file based on the input sequences.
the output files are saved in the specified output directory.
Usage:
   python full_length_make_gtf_ref.py -i <input_fasta> -o <output_directory>
Arguments:
   -i, --input     Input FASTA file or directory containing FASTA files (required)
   -o, --output    Output directory path (required)
Example:
   python full_length_make_gtf_ref.py -i input.fasta -o output/
Author: Shengyao Zhang
Date: 2024-12-19
Version: 1.0
"""
import argparse
import os
import sys

class BaseParser:
    def __init__(self, input_fp="", output_dir="") -> None:
        self.input_fp = input_fp
        self.output_dir = output_dir
        
        self._check_input_path()
        self._create_output_directory()
    
    def _check_input_path(self):
        """
        Check if the input path exists. If it's a directory, look for .fa or .fasta files.
        """
        if not os.path.exists(self.input_fp):
            raise FileNotFoundError(f"Input path does not exist: {self.input_fp}")
        
        if os.path.isdir(self.input_fp):
            fasta_files = [f for f in os.listdir(self.input_fp) if f.endswith(('.fa', '.fasta'))]
            if not fasta_files:
                raise FileNotFoundError(f"No FASTA file found in directory: {self.input_fp}")
            self.input_fp = os.path.join(self.input_fp, fasta_files[0])
            print(f"Using FASTA file: {self.input_fp}")
        
        if not os.access(self.input_fp, os.R_OK):
            raise PermissionError(f"Cannot read input file: {self.input_fp}")
    
    def _create_output_directory(self):
        """
        Create the output directory. If the directory already exists, skip creation.
        """
        try:
            os.makedirs(self.output_dir, exist_ok=True)
            print(f"Output directory: {self.output_dir}")
        except Exception as e:
            raise Exception(f"Failed to create output directory: {str(e)}")
    
    def read_file(self):
        raise NotImplementedError("Subclass must implement read_file method")
    
    def write_sequence(self):
        raise NotImplementedError("Subclass must implement write_sequence method")
    
    def run(self):
        self.read_file()

class FtoGnGtf(BaseParser):
    def __init__(self, input_fp="", output_dir="") -> None:
        super().__init__(input_fp, output_dir)
        self.current_sequence_id = None
        self.current_sequence = []
        self.genome_fasta_fp = os.path.join(self.output_dir, "output.fasta")
        self.gtf_fp = os.path.join(self.output_dir, "output.gtf")
        
        self._clear_existing_output_files()

    def _clear_existing_output_files(self):
        """
        Clear existing output files.
        """
        for fp in [self.genome_fasta_fp, self.gtf_fp]:
            if os.path.exists(fp):
                os.remove(fp)
                print(f"Removed existing file: {fp}")

    def read_file(self):
        try:
            with open(self.input_fp, "r") as fts:
                for line_num, line in enumerate(fts, 1):
                    try:
                        if line.startswith('>'):
                            if self.current_sequence_id is not None:
                                self.write_sequence()
                            self.current_sequence_id = line[1:].split()[0].strip()
                            self.current_sequence = []
                        else:
                            self.current_sequence.append(line.strip())
                    except Exception as e:
                        print(f"Error processing line {line_num}: {str(e)}")
                        continue
                
                if self.current_sequence_id is not None:
                    self.write_sequence()
        except Exception as e:
            raise Exception(f"Error reading file: {str(e)}")

    def write_sequence(self):
        try:
            sequence = ''.join(self.current_sequence)
            attrib_record = {
                "gene_id": self.current_sequence_id,
                "transcript_id": self.current_sequence_id,
                "exon_number": "1",
                "gene_name": self.current_sequence_id,
                "gene_biotype": "protein_coding"
            }
            attr_rec = [f'{k} "{v}"' for k, v in attrib_record.items()]
            attr_rec = "; ".join(attr_rec)
            start_p = "1"
            end_p = str(len(sequence))

            record = [
                self.current_sequence_id,
                ".",
                "cds",
                start_p,
                end_p,
                ".",
                "+",
                ".",
                attr_rec
            ]
            record_str = "\t".join(record)

            with open(self.genome_fasta_fp, "a") as fasta_fp:
                fasta_fp.write(f">{self.current_sequence_id}\n{sequence}\n")

            with open(self.gtf_fp, "a") as gtf_fp:
                gtf_fp.write(f"{record_str}\n")
                
        except Exception as e:
            print(f"Error writing sequence ({self.current_sequence_id}): {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='Convert FASTA to GTF and new FASTA.')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file or directory containing FASTA files')
    parser.add_argument('-o', '--output', required=True, help='Output directory path')
    args = parser.parse_args()

    try:
        f_to_gngtf = FtoGnGtf(args.input, args.output)
        f_to_gngtf.run()
        print("Conversion completed!")
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
