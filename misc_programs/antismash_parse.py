from Bio import SeqIO
from annot_fungi import *
import os, shutil, glob
import pandas as pd

def get_by_end(path, end):
    for file in os.listdir(path):
        if file.endswith(end):
            return f"{path}/{file}"

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

def main():
    target='tw_fun_assem'
    for path in os.listdir(target):
        FUNG = fungi(f"{target}/{path}")
        annotate_result_path = f"{FUNG.annot_dir}"
        antismash_path = f"{get_by_end(annotate_result_path, 'antismash_folder')}/"
        gbk_path=find_largest_gbk_file(antismash_path)
        gbk_folder='antismash_gbk'
        os.makedirs(gbk_folder, exist_ok=True)

        shutil.copy(gbk_path, f"{gbk_folder}/{FUNG.basename}.gbk")


    for files in glob.glob(f"{gbk_folder}/*.gbk"):
        tsv_folder= "tsv_files"
        os.makedirs(tsv_folder, exist_ok=True)
        out_files = files.replace(".gbk","_output.tsv")
        cluster_out = open(out_files, "w")


    # Extract Cluster info, write to file
        for seq_record in SeqIO.parse(files, "genbank"):
            for seq_feat in seq_record.features:
                if seq_feat.type == "protocluster":
                    cluster_number = seq_feat.qualifiers["protocluster_number"][0].replace(" ","_").replace(":","")
                    cluster_type = seq_feat.qualifiers["product"][0]
                    cluster_out.write("#"+cluster_number+"\tCluster Type:"+cluster_type+"\n")

    df = {}

    for fileName in glob.glob(f"{gbk_folder}/*.tsv"):  # Prefix "*.tsv" with the target directory
        genomeIdentifier = fileName.split("/")[-1].replace(".tsv", "")
        for line in open(fileName):
            line = line.strip()
            line = line.split("\t")[-1].replace("Cluster Type:", "")
            if (not line in df):
                df[line] = {}
            df[line][genomeIdentifier] = df[line].get(genomeIdentifier, 0) + 1
    df = pd.DataFrame(df)
    df = df.fillna(0)  # To handle missing value
    df.to_csv("table.tsv", sep="\t")

if __name__ == "__main__":
    main()