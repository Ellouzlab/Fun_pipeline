import argparse, os, shutil
import subprocess as sp
from multiprocessing import cpu_count

def arg_parser():
    parser = argparse.ArgumentParser(description="Run funannotate pipeline")
    parser.add_argument("-i", "--input", help="input fasta", required=True)
    parser.add_argument("-o", "--output", help="Output directory", required=True)
    parser.add_argument("-s", "--species", help="Species name in the form 'Genus species''", required=True)
    parser.add_argument("-t", "--threads", help="Number of threads (default: Max threads - 1)", required=False, default=cpu_count()-1)
    return parser.parse_args()


def clean(input_path: str):
    output_path=input_path.replace(".fasta", "_clean.fasta")
    command = f"funannotate clean -i {input_path} -o {output_path}"
    if not os.path.exists(output_path):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

def sort(input_path: str):
    output_path=input_path.replace("_clean.fasta", "_sort.fasta")
    command = f"funannotate sort -i {input_path} -o {output_path}"
    if not os.path.exists(output_path):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

def mask(input_path: str, species: str, threads: int):
    output_path=input_path.replace("_sort.fasta", "_mask.fasta")
    command = f"funannotate mask -i {input_path} -o {output_path} -s '{species}' --cpus {threads}"
    if not os.path.exists(output_path):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

def predict(input_path: str, species: str, threads: int):
    output_path=input_path.replace("_mask.fasta", "_predict_folder")
    command = f"funannotate predict -i {input_path} -o {output_path} -s '{species}' --cpus {threads}"
    if not os.path.exists(output_path):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

def iprscan(input_path: str, threads: int):
    output_path=input_path.replace("_predict_folder", "_ipr.xml")
    command = f"funannotate iprscan -i {input_path} -o {output_path} -m docker --cpus {threads}"
    if not os.path.exists(output_path):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

def antismash(input_path: str):
    result_dir_path=f"{input_path}/predict_results"
    gff_path="unknown"
    gbk_path="unknown"
    output_path=input_path.replace("_predict_folder", "_antismash_folder")

    for file in os.listdir(result_dir_path):
        if file.endswith(".gbk"):
            gbk_path=os.path.join(result_dir_path, file)
        if file.endswith(".gff3"):
            gff_path=os.path.join(result_dir_path, file)
    if gff_path=="unknown" or gbk_path=="unknown":
        raise Exception("Antismash failed to run because predict results gbk or gff were not found")
    
    with open("tmp_change_env.sh", "w") as file:
        file.write(f"conda activate antismash\n")
        file.write(f"antismash --taxon fungi --genefinding-gff3 {gff_path} --output-dir {output_path} {gbk_path}\n")
        file.write(f"conda activate funannotate\n")
    
    if not os.path.exists(output_path):
        print("Running antismash")
        sp.run("bash -i tmp_change_env.sh", shell=True, executable="/bin/bash")
    
    os.remove("tmp_change_env.sh")
    return output_path

def signalp(input_path: str):
    result_dir_path=f"{input_path}/predict_results"
    prot_path="unknown"
    output_path=input_path.replace("_predict_folder", "_signalp_folder")

    for file in os.listdir(result_dir_path):
        if file.endswith(".proteins.fa"):
            prot_path=os.path.join(result_dir_path, file)
    if prot_path=="unknown":
        raise Exception("signalp failed to run because predict results proteins.fa was not found")
    command=f"signalp6 --fastafile {prot_path} --organism euk --mode fast --output_dir {output_path}"
    print(prot_path)
    if not os.path.exists(output_path):
        sp.run(command, shell=True, executable="/bin/bash")

    return output_path

def find_largest_gbk_file(directory):
    largest_file = None
    largest_size = 0

    for filename in os.listdir(directory):
        if filename.endswith(".gbk"):
            filepath = os.path.join(directory, filename)
            size = os.path.getsize(filepath)

            if size > largest_size:
                largest_size = size
                largest_file = filepath

    return largest_file

def annotate(input_path: str, species: str, iprscan_path: str, antismash_dir: str, signalp_dir:str, threads: int):
    output_path=input_path.replace("_predict_folder", "_annotate_folder")

    antismash_gbk=find_largest_gbk_file(antismash_dir)
    signalp_gff3=f"{signalp_dir}/prediction_results.txt"

    command = f"funannotate annotate -i {input_path} -o {output_path} --cpus {threads} -s '{species}'"
    command += f" --antismash {antismash_gbk} --signalp {signalp_gff3} --iprscan {iprscan_path}"
    if not os.path.exists(output_path):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

def main():
    args=arg_parser()
    os.makedirs(args.output, exist_ok=True)
    new_input_path=f"{args.output}/{args.input.split('/')[-1]}"
    shutil.copy(args.input, new_input_path)
    
    #Run Pipeline
    clean_path=clean(new_input_path)
    sort_path=sort(clean_path)
    mask_path=mask(sort_path, args.species, args.threads)
    predict_path=predict(mask_path, args.species, args.threads)
    ipr_path=iprscan(predict_path, args.threads)
    antismash_path=antismash(predict_path)
    signalp_path=signalp(predict_path)
    annotate(predict_path, args.species, ipr_path, antismash_path, signalp_path, args.threads)

if __name__ == "__main__":
    main()