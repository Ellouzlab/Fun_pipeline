import os, shutil
def main():
    annotations="annotations"
    for fungidir in os.listdir(annotations):
        fungidir_path=f"{annotations}/{fungidir}"

        base='_'.join(fungidir.split("_")[:-1])
        species=fungidir.split("_")[1]
        genus = fungidir.split("_")[0]

        annotdir_path=f"{fungidir_path}/{base}_predict_folder/annotate_results"

        fasta_path=f"{annotdir_path}/{species}.scaffolds.fa"
        alt_fasta_path=f"{annotdir_path}/{genus}_{species}.scaffolds.fa"
        gff3_path=f"{annotdir_path}/{species}.gff3"
        tbl_path=f"{annotdir_path}/{species}.tbl"
        alt_tbl_path=f"{annotdir_path}/{genus}_{species}.tbl"
        alt_gff3_path=f"{annotdir_path}/{genus}_{species}.gff3"
        gbk_path = f"{annotdir_path}/{species}.gbk"
        alt_gbk_path = f"{annotdir_path}/{genus}_{species}.gbk"

        outdir=f"annotation_folders/{base}"
        os.makedirs(outdir,exist_ok=True)

        try:
            shutil.copy(fasta_path,f"{outdir}/{base}.fsa")
            shutil.copy(gff3_path,f"{outdir}/{base}.gff3")
            shutil.copy(gbk_path,f"{outdir}/{base}.gbk")
            shutil.copy(tbl_path,f"{outdir}/{base}.tbl")
        except:
            shutil.copy(alt_fasta_path,f"{outdir}/{base}.fsa")
            shutil.copy(alt_gff3_path,f"{outdir}/{base}.gff3")
            shutil.copy(alt_gbk_path,f"{outdir}/{base}.gbk")
            shutil.copy(alt_tbl_path,f"{outdir}/{base}.tbl")

if __name__ == '__main__':
    main()