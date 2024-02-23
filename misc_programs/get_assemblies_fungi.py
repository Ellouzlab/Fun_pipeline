import os, shutil
from annot_fungi import *

def arg_parse():
    parser = argparse.ArgumentParser(description="Get assemblies")
    parser.add_argument("-i", "--indir", help="input dir", required=True)
    parser.add_argument('-e', '--extension', help='extension of file', required=True)
    parser.add_argument('-o', '--outdir', help='output dir', required=True)
    parser.add_argument('-n', '--new_extension', help='New Extension for files', required=False, default="")
    return parser.parse_args()

def main():
    args = arg_parse()
    os.makedirs(args.outdir, exist_ok=True)
    for path in os.listdir(args.indir):
        FUNG = fungi(f"{args.indir}/{path}")
        annotate_result_path = f"{FUNG.predict_dir}/annotate_results"


        for file in os.listdir(annotate_result_path):
            if file.endswith(args.extension):
                file_path= f"{annotate_result_path}/{file}"
                print(file_path)
                new_file_path=f"{args.outdir}/{FUNG.basename}.{args.new_extension}"
                shutil.copy(file_path, new_file_path)

if __name__ == "__main__":
    main()