import argparse, os, shutil
import subprocess as sp
from multiprocessing import cpu_count
import pandas as pd
from datetime import datetime
'''from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF'''


def arg_parser():
    parser = argparse.ArgumentParser(description="Run funannotate pipeline")
    parser.add_argument("-i", "--input", help="input fasta", required=True)
    parser.add_argument("-o", "--output", help="Output directory", required=True)
    parser.add_argument("-s", "--species", help="Species name in the form 'Genus species''", required=True)
    parser.add_argument("-t", "--threads", help="Number of threads (default: Max threads - 1)", required=False, default=cpu_count()-1)
    parser.add_argument("--export", help="Path to fun_exports.sh file with paths", required=False, default="/home/sulman/oualid_fungi/funannotate_pipeline/fun_exports.sh")
    return parser.parse_args()

def running_message(function):
    def wrapper(*args, **kwargs):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        args_repr = [repr(a) for a in args]
        kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
        signature = ", ".join(args_repr + kwargs_repr)
        print("----------------------------------------")
        print(f"Time: {current_time} - Running {function.__name__}\n\nUsing inputs:\n{function.__name__}({signature})\n")
        result = function(*args, **kwargs)

        now2 = datetime.now()
        current_time2 = now2.strftime("%H:%M:%S")
        print(f"\nTime: {current_time2} - {function.__name__} Completed\n\n")
        return result
    return wrapper

def export_path(export_file_path: str):
    print(f"\n\nExporting paths to programs and databases. Please remember to check {export_file_path}\n")
    print('Paths are the following:')
    command = f"source {export_file_path}"
    sp.run(command, shell=True, executable="/bin/bash")

    with open(export_file_path) as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("export"):
            new_line=line.replace('  ', ' ')
            new_line=new_line.replace("export ", "echo $")
            new_line=new_line.split('=')[0]
            print(f"\n{new_line.split(' ')[1]}:")
            sp.run(new_line, shell=True, executable="/bin/bash")

@running_message
def clean(input_path: str):
    output_path=input_path.replace(".fasta", "_clean.fasta")
    command = f"funannotate clean -i {input_path} -o {output_path}"
    if not os.path.exists(output_path):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

@running_message
def sort(input_path: str):
    output_path=input_path.replace("_clean.fasta", "_sort.fasta")
    command = f"funannotate sort -i {input_path} -o {output_path}"
    if not os.path.exists(output_path):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

@running_message
def mask(input_path: str, species: str, threads: int):
    output_path=input_path.replace("_sort.fasta", "_mask.fasta")
    command = f"funannotate mask -i {input_path} -o {output_path} -s '{species}' --cpus {threads}"
    if not os.path.exists(output_path):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

@running_message
def predict(input_path: str, species: str, threads: int):
    output_path=input_path.replace("_mask.fasta", "_predict_folder")
    command = f"funannotate predict -i {input_path} -o {output_path} -s '{species}' --cpus {threads}"
    if not os.path.exists(output_path):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

@running_message
def iprscan(input_path: str, threads: int):
    output_path=input_path.replace("_predict_folder", "_ipr.xml")
    command = f"funannotate iprscan -i {input_path} -o {output_path} -m docker --cpus {threads}"
    if not os.path.exists(output_path):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

@running_message
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

@running_message
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

@running_message
def annotate(input_path: str, species: str, iprscan_path: str, antismash_dir: str, signalp_dir:str, threads: int):
    output_path=input_path.replace("_predict_folder", "_annotate_folder")

    antismash_gbk=find_largest_gbk_file(antismash_dir)
    signalp_gff3=f"{signalp_dir}/prediction_results.txt"

    command = f"funannotate annotate -i {input_path} -o {output_path} --cpus {threads} -s '{species}'"
    command += f" --antismash {antismash_gbk} --signalp {signalp_gff3} --iprscan {iprscan_path}"
    if not os.path.exists(f"{input_path}/annotate_results"):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

@running_message
def barrnap(input_path: str, predict_dir, threads: int):
    output_path=predict_dir.replace("_predict_folder", "_barrnap_folder")
    os.makedirs(output_path, exist_ok=True)
    command = f"barrnap --kingdom euk --threads {threads} {input_path} > {output_path}/barrnap.gff"
    if not os.path.exists(f"{output_path}/barrnap.gff"):
        sp.run(command, shell=True, executable="/bin/bash")
    return output_path

def modified_barrnap(input_path):
    barrnap_path=f"{input_path}/barrnap.gff"
    modified_barrnap_path=f"{input_path}/modified_barrnap.gff3"
    if not os.path.exists(modified_barrnap_path):
        barrnap_df=pd.read_csv(barrnap_path, sep="\t", skiprows=1)
        barrnap_df.iloc[:,8]=barrnap_df.iloc[:,8].str.split("note").str[0]
        barrnap_df.iloc[:,5]='.'
        barrnap_df.to_csv(modified_barrnap_path, sep="\t", index=False, header=False)
        with open(modified_barrnap_path, "r") as file:
            data=file.read()
        new_data=f"##gff-version 3\n{data}"
        with open(modified_barrnap_path, "w") as file:
            file.write(new_data)


def add_gff_to_gbk(gff_file, gbk_file, output_file):
    # Read the GenBank file
    gbk_records = SeqIO.to_dict(SeqIO.parse(gbk_file, "genbank"))

    # Initialize a dictionary to store GFF information
    gff_data = {}

    # Parse the GFF file and store the relevant information
    with open(gff_file, "r") as gff_handle:
        for rec in GFF.parse(gff_handle):
            for feature in rec.features:
                if feature.type == "rRNA":
                    locus_tag = feature.qualifiers["Name"][0]
                    product = feature.qualifiers["product"][0]
                    start = feature.location.start.position
                    end = feature.location.end.position
                    strand = "+" if feature.location.strand == 1 else "-"

                    if rec.id not in gff_data:
                        gff_data[rec.id] = []

                    gff_data[rec.id].append({
                        "locus_tag": locus_tag,
                        "product": product,
                        "start": start,
                        "end": end,
                        "strand": strand
                    })

    # Update the GenBank records with GFF information
    for record_id, gff_features in gff_data.items():
        if record_id in gbk_records:
            for feature_info in gff_features:
                gbk_record = gbk_records[record_id]
                new_feature = SeqFeature()  # Use SeqFeature from Bio.SeqFeature
                new_feature.location = FeatureLocation(  # Use FeatureLocation from Bio.SeqFeature
                    start=feature_info["start"],
                    end=feature_info["end"],
                    strand=1 if feature_info["strand"] == "+" else -1
                )
                new_feature.type = "rRNA"
                new_feature.qualifiers["locus_tag"] = feature_info["locus_tag"]
                new_feature.qualifiers["product"] = feature_info["product"]

                gbk_record.features.append(new_feature)

    # Write the updated GenBank records to a new file
    with open(output_file, "w") as output_handle:
        for record_id, gbk_record in gbk_records.items():
            SeqIO.write(gbk_record, output_handle, "genbank")

def extract_sequences_and_tbl_from_gbk(gbk_file, output_prefix):
    records = SeqIO.parse(gbk_file, "genbank")

    ffn_filename = f"{output_prefix}.ffn"
    fna_filename = f"{output_prefix}.fna"
    faa_filename = f"{output_prefix}.faa"
    tbl_filename = f"{output_prefix}.tbl"

    with open(ffn_filename, "w") as ffn_file, \
         open(fna_filename, "w") as fna_file, \
         open(faa_filename, "w") as faa_file, \
         open(tbl_filename, "w") as tbl_file:

        for record in records:
            for feature in record.features:
                if feature.type == "CDS":
                    locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                    protein_id = feature.qualifiers.get("protein_id", ["unknown"])[0]
                    start = feature.location.start.position + 1  # Convert 0-based to 1-based
                    end = feature.location.end.position
                    strand = "+" if feature.location.strand == 1 else "-"

                    ffn_file.write(f">{locus_tag}\n{feature.location.extract(record).seq}\n")
                    fna_file.write(f">{locus_tag}\n{record.seq}\n")

                    if "translation" in feature.qualifiers:
                        faa_file.write(f">{protein_id}\n{feature.qualifiers['translation'][0]}\n")

                    # Write information to .tbl file
                    tbl_file.write(f"{start}\t{end}\tCDS\n")
                    tbl_file.write(f"\t\t\tlocus_tag\t{locus_tag}\n")
                    tbl_file.write(f"\t\t\tprotein_id\t{protein_id}\n")
                    tbl_file.write(f"\t\t\tproduct\t{feature.qualifiers['product'][0]}\n")
                    tbl_file.write(f"\t\t\tprotein_id\t{protein_id}\n")
                    tbl_file.write(f"\t\t\ttransl_table\t11\n")
                    tbl_file.write(f"\t\t\tproduct\t{feature.qualifiers['product'][0]}\n")

@running_message
def finalizing_annotations(barrnap_folder, predict_path, output_path):
    results_folder=f"{output_path}/results"
    os.makedirs(results_folder, exist_ok=True)
    print("Finalizing Results...")

    gff_path=f"{barrnap_folder}/modified_barrnap.gff3"

    annotation_path=f"{predict_path}/annotate_results"
    for file in os.listdir(annotation_path):
        if file.endswith(".gbk"):
            gbk_file=file

    gbk_path=f"{annotation_path}/{gbk_file}"

    output_gbk=f"{results_folder}/annotations.gbk"
    add_gff_to_gbk(gff_path, gbk_path, output_gbk)
    extract_sequences_and_tbl_from_gbk(output_gbk, f"{results_folder}/annotations")


def main():
    args=arg_parser()
    os.makedirs(args.output, exist_ok=True)
    new_input_path=f"{args.output}/{args.input.split('/')[-1]}"
    shutil.copy(args.input, new_input_path)
    
    #Run Main Pipeline
    export_path(args.export)
    clean_path=clean(new_input_path)
    sort_path=sort(clean_path)
    mask_path=mask(sort_path, args.species, args.threads)
    predict_path=predict(mask_path, args.species, args.threads)
    ipr_path=iprscan(predict_path, args.threads)
    antismash_path=antismash(predict_path)
    signalp_path=signalp(predict_path)
    annotate(predict_path, args.species, ipr_path, antismash_path, signalp_path, args.threads)
    '''
    #Run barrnap
    annot_folder=f"{predict_path}/annotate_results"
    for file in os.listdir(annot_folder):
        if file.endswith(".scaffolds.fa"):
            fasta_path=f"{annot_folder}/{file}"

    barrnap_folder=barrnap(fasta_path, predict_path, args.threads)
    modified_barrnap(barrnap_folder)
    finalizing_annotations(barrnap_folder, predict_path, args.output)'''




if __name__ == "__main__":
    main()