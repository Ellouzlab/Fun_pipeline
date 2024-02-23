import os
import subprocess as sp
from BCBio import GFF

def fix_gbk(gbk_path):
    with open(gbk_path, "r") as f:
        lines=f.readlines()
    to_remove = []
    filtered_lines=[line for line in lines if not "EC_number" in line]
    gbk_txt=''.join(filtered_lines)
    with open(gbk_path, "w") as f:
        f.write(gbk_txt)
def fix_tbl(tbl_path):
    with open(tbl_path, "r") as f:
        lines=f.readlines()

    gene_lines=[]
    for i in range(len(lines)):
        try:
            if lines[i].split()[2]=="gene":
                gene_lines.append(i)
        except:
            pass
    hypo_dict = {}
    to_remove = []
    for k in range(len(gene_lines)-1):
        start=gene_lines[k]
        end=gene_lines[k+1]

        for index in range(start,end):
            if "hypothetical" in lines[index]:
                hypo_dict[index]=[start, end]
                break


        for key in hypo_dict:
            start=hypo_dict[key][0]
            end=hypo_dict[key][1]
            for index in range(start,end):
                if "EC_number" in lines[index]:
                    if not index in to_remove:
                        to_remove.append(index)

    filtered_lines=[lines[i] for i in range(len(lines)) if not i in to_remove]
    tbl_txt=''.join(filtered_lines)

    with open(tbl_path, "w") as f:
        f.write(tbl_txt)


def get_locustag(gff3_path):
    in_handle = open(gff3_path)
    for record in GFF.parse(in_handle):
        for feature in record.features:
            if feature.type=="gene":
                return feature.id.split('_')[0]

def run_tbl2asn(species, isolate, fungidir, gff_path, tbl_path, locustag):
    cmd=f"table2asn.linux64 -M n -J -c w -euk -t template.sbt -gaps-min 10 -l paired-ends -j '[organism={species}] [isolate={isolate}]' "
    cmd+=f"-indir {fungidir} -Z -locus-tag-prefix {locustag} -outdir {fungidir} -f {tbl_path}"
    sp.run(cmd,shell=True)

def main():
    annotations="annotation_folders"
    for fungidir in os.listdir(annotations):
        fungidir_path=f"{annotations}/{fungidir}"
        species=" ".join(fungidir.split("_")[0:2])
        isolate=fungidir.split("_")[-1]
        print(species, isolate)

        tbl_path=f"{fungidir_path}/{fungidir}.tbl"
        gff3_path=f"{fungidir_path}/{fungidir}.gff3"
        gbk_path=f"{fungidir_path}/{fungidir}.gbk"
        fix_gbk(tbl_path)
        locustag=get_locustag(gff3_path)
        #fix_tbl(tbl_path)

        run_tbl2asn(species, isolate, fungidir_path, gff3_path, tbl_path, locustag)


if __name__ == '__main__':
    main()