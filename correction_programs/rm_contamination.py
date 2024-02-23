from Bio import SeqIO
from BCBio import GFF

def read_txt(file):
    with open(file, "r") as f:
        return f.read()

def write_txt(file, text):
    with open(file, "w") as f:
        f.write(text)
def read_fasta(file):
    return SeqIO.parse(file, "fasta")


def write_fasta(file, records):
    with open(file, "w") as f:
        SeqIO.write(records, f, "fasta")
def filter_and_write_gff3(input_file, output_file, ids_to_remove):
    lines=read_txt(input_file).split("\n")
    n_lines=[]
    for line in lines:
        found=0
        for id in ids_to_remove:
            if id in line:
                found=1
        if found==1:
            continue
        else:
            n_lines.append(line)
    write_txt(output_file,'\n'.join(n_lines))


def main():
    lines=read_txt("diatrype_contamination.txt").split("\n")

    scaffolds=[]
    for line in lines:
        scaff=line.split('\t')[0]
        if not scaff=='':
            scaffolds.append(scaff)

    fasta_path="annotation_folders/Diatrype_stigma_strain_M11_M66-122/Diatrype_stigma_strain_M11_M66-122.fsa"
    records=read_fasta(fasta_path)
    filtered_records=[record for record in records if record.id not in scaffolds]
    write_fasta("annotation_folders/Diatrype_stigma_strain_M11_M66-122/Diatrype_stigma_strain_M11_M66-122.fsa", filtered_records)

    gff3_path="annotation_folders/Diatrype_stigma_strain_M11_M66-122/Diatrype_stigma_strain_M11_M66-122.gff3"
    mod_gff3_path="annotation_folders/Diatrype_stigma_strain_M11_M66-122/Diatrype_stigma_strain_M11_M66-122.gff3"

    filter_and_write_gff3(gff3_path, mod_gff3_path, scaffolds)


if __name__ == '__main__':
    main()