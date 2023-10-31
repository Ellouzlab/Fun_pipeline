import argparse, os
import subprocess as sp
from multiprocessing import cpu_count
def arg_parser():
    parser = argparse.ArgumentParser(description="Run funannotate pipeline")
    parser.add_argument("-i", "--input", help="input dir", required=True)
    parser.add_argument("-o", "--output", help="Output directory", required=True)
    parser.add_argument("-t", "--threads", help="Number of threads (default: Max threads - 1)", required=False, default=cpu_count()-1)
    return parser.parse_args()

def main():
    args=arg_parser()
    os.makedirs(args.output, exist_ok=True)

    for file in os.listdir(args.input):
        file_path=os.path.join(args.input, file)
        if file.endswith(".fasta"):
            output_annotation_path=os.path.join(args.output, file.replace(".fasta", "_annotation"))
            
            if not os.path.exists(output_annotation_path):
                os.makedirs(output_annotation_path)
                species=' '.join(file.split("_")[1:2])
                command = f"python funannotate_pipeline/pipeline.py -i {file_path} -o {output_annotation_path} -t {args.threads} -s '{species}'"

                sp.run(command, shell=True, executable="/bin/bash")


if __name__ == "__main__":
    main()