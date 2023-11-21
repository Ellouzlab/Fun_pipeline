import argparse, os, sys,  shutil
from Bio import SeqIO


def arg_parser():
    parser = argparse.ArgumentParser(description="Rename fasta files downloaded from NCBI genbank for the purpose of annotation")
    parser.add_argument("-i", "--indir", help="Dir with other dir fasta", required=True)
    parser.add_argument("-o", "--output", help="Output directory", required=True)
    return parser.parse_args()

def fasta_reader(file):
    for record in SeqIO.parse(file, "fasta"):
        return (f"{'_'.join(record.description.split(' ')[1:5])}.fasta").replace("/", "_")

def main():
    args = arg_parser()
    os.makedirs(args.output, exist_ok=True)
    for stuff in os.listdir(args.indir):
        stuff_path=os.path.join(args.indir, stuff)
        if os.path.isdir(stuff_path):
            for file in os.listdir(stuff_path):
                if file.endswith(".fna"):
                    file_path=os.path.join(stuff_path, file)
                    new_name=fasta_reader(file_path)
                    shutil.copy(file_path, os.path.join(args.output, new_name))
if __name__ == "__main__":
    main()