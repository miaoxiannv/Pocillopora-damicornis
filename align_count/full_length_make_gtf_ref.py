import argparse
import os
import sys

class BaseParser:
    def __init__(self, input_fp="", output_dir="") -> None:
        self.input_fp = input_fp
        self.output_dir = output_dir
        
        # 检查输入路径
        if not os.path.exists(input_fp):
            raise FileNotFoundError(f"输入路径不存在: {input_fp}")
        
        # 如果输入是目录，查找.fa或.fasta文件
        if os.path.isdir(input_fp):
            fasta_files = [f for f in os.listdir(input_fp) if f.endswith(('.fa', '.fasta'))]
            if not fasta_files:
                raise FileNotFoundError(f"在目录中未找到FASTA文件: {input_fp}")
            self.input_fp = os.path.join(input_fp, fasta_files[0])
            print(f"使用FASTA文件: {self.input_fp}")
        
        # 检查输入文件权限
        if not os.access(self.input_fp, os.R_OK):
            raise PermissionError(f"无法读取输入文件: {self.input_fp}")
        
        # 创建输出目录
        try:
            os.makedirs(output_dir, exist_ok=True)
            print(f"输出目录: {output_dir}")
        except Exception as e:
            raise Exception(f"创建输出目录失败: {str(e)}")
    
    def read_file(self):
        raise NotImplementedError("子类必须实现read_file方法")
    
    def write_sequence(self):
        raise NotImplementedError("子类必须实现write_sequence方法")
    
    def run(self):
        self.read_file()

class FtoGnGtf(BaseParser):
    def __init__(self, input_fp="", output_dir="") -> None:
        super().__init__(input_fp, output_dir)
        self.current_sequence_id = None
        self.current_sequence = []
        self.genome_fasta_fp = os.path.join(self.output_dir, "output.fasta")
        self.gtf_fp = os.path.join(self.output_dir, "output.gtf")
        
        # 清除现有输出文件
        for fp in [self.genome_fasta_fp, self.gtf_fp]:
            if os.path.exists(fp):
                os.remove(fp)
                print(f"已删除现有文件: {fp}")

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
                        print(f"处理第{line_num}行时出错: {str(e)}")
                        continue
                
                if self.current_sequence_id is not None:
                    self.write_sequence()
        except Exception as e:
            raise Exception(f"读取文件时出错: {str(e)}")

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
            print(f"写入序列时出错 ({self.current_sequence_id}): {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='Convert FASTA to GTF and new FASTA.')
    parser.add_argument('-i', '--input', required=True, help='输入FASTA文件或包含FASTA文件的目录路径')
    parser.add_argument('-o', '--output', required=True, help='输出文件夹路径')
    args = parser.parse_args()

    try:
        f_to_gngtf = FtoGnGtf(args.input, args.output)
        f_to_gngtf.run()
        print("转换完成!")
    except Exception as e:
        print(f"错误: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
