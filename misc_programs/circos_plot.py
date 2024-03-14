from pycirclize import Circos
from pycirclize.parser import Genbank
from pycirclize.utils import ColorCycler
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import threading, argparse, shutil, os
from multiprocessing import cpu_count
import pandas as pd
import concurrent.futures
from Bio.Blast.Applications import NcbiblastnCommandline


def argparser():
    parser = argparse.ArgumentParser(description='Circos plot')
    parser.add_argument('-g', '--genbank', type=str, help='Genbank file', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file')
    parser.add_argument('-t', '--threads', type=int, required=False, default=cpu_count()-1, help='Number of threads to use')
    return parser.parse_args()

def calc_genome_gc_content(seq) -> float:
    try:
        gc_content = SeqUtils.gc_fraction(seq) * 100
    except AttributeError:
        gc_content = SeqUtils.GC(seq)
    return gc_content

def calc_gc_content(seq, window_size=10000, step_size=5000):
    #stolen from moshi4
    pos_list, gc_content_list = [], []
    pos_list = list(range(0, len(seq), step_size)) + [len(seq)]
    for pos in pos_list:
        window_start_pos = pos - int(window_size / 2)
        window_end_pos = pos + int(window_size / 2)
        window_start_pos = 0 if window_start_pos < 0 else window_start_pos
        window_end_pos = len(seq) if window_end_pos > len(seq) else window_end_pos

        subseq = seq[window_start_pos:window_end_pos]
        try:
            gc_content = SeqUtils.gc_fraction(subseq) * 100
        except AttributeError:
            # For backward compatibility, biopython <= 1.79
            gc_content = SeqUtils.GC(subseq)  # type: ignore
        gc_content_list.append(gc_content)
    return (np.array(pos_list), np.array(gc_content_list))

def calc_gc_skew(seq, window_size=10000, step_size=5000) -> tuple[np.ndarray, np.ndarray]:
    pos_list, gc_skew_list = [], []
    pos_list = list(range(0, len(seq), step_size)) + [len(seq)]
    for pos in pos_list:
        window_start_pos = pos - int(window_size / 2)
        window_end_pos = pos + int(window_size / 2)
        window_start_pos = 0 if window_start_pos < 0 else window_start_pos
        window_end_pos = len(seq) if window_end_pos > len(seq) else window_end_pos

        subseq = seq[window_start_pos:window_end_pos]
        g = subseq.count("G") + subseq.count("g")
        c = subseq.count("C") + subseq.count("c")
        try:
            skew = (g - c) / float(g + c)
        except ZeroDivisionError:
            skew = 0.0
        gc_skew_list.append(skew)

    return (np.array(pos_list), np.array(gc_skew_list))

def count_special_sequences(seq, window_size=10000, step_size=5000) -> tuple[np.ndarray, np.ndarray]:
    pos_list, sequence_count_list = [], []
    pos_list = list(range(0, len(seq), step_size)) + [len(seq)]
    for pos in pos_list:
        window_start_pos = pos - int(window_size / 2)
        window_end_pos = pos + int(window_size / 2)
        window_start_pos = 0 if window_start_pos < 0 else window_start_pos
        window_end_pos = len(seq) if window_end_pos > len(seq) else window_end_pos

        subseq = seq[window_start_pos:window_end_pos].upper()  # Convert to uppercase to avoid case-sensitivity issues
        count_ttaggg = subseq.count("TTAGGG")
        count_ccctaa = subseq.count("CCCTAA")
        total_count = count_ttaggg + count_ccctaa
        sequence_count_list.append(total_count)

    return (np.array(pos_list), np.array(sequence_count_list))

def add_gc_content_skew_tracks(circos, seqid2seq):
    for sector in circos.sectors:
        seq = seqid2seq[sector.name]

        if len(seq) > 100000:
            pos_list, gc_contents = calc_gc_content(seq)
            gc_content_track = sector.add_track((63, 73), r_pad_ratio=0.1)
            gc_skew_track = sector.add_track((53, 63), r_pad_ratio=0.1)
            telomere_track = sector.add_track((43, 53), r_pad_ratio=0.1)

            #GC content
            gc_contents = gc_contents - calc_genome_gc_content(seq)
            positive_gc_contents = np.where(gc_contents > 0, gc_contents, 0)
            negative_gc_contents = np.where(gc_contents < 0, gc_contents, 0)
            abs_max_gc_content = np.max(np.abs(gc_contents))
            vmin, vmax = -abs_max_gc_content, abs_max_gc_content
            gc_content_track.fill_between(
                pos_list, positive_gc_contents, 0, vmin=vmin, vmax=vmax, color="black"
            )
            gc_content_track.fill_between(
                pos_list, negative_gc_contents, 0, vmin=vmin, vmax=vmax, color="grey"
            )

            pos_list, gc_skews = calc_gc_skew(seq)
            positive_gc_skews = np.where(gc_skews > 0, gc_skews, 0)
            negative_gc_skews = np.where(gc_skews < 0, gc_skews, 0)
            abs_max_gc_skew = np.max(np.abs(gc_skews))
            vmin, vmax = -abs_max_gc_skew, abs_max_gc_skew
            gc_skew_track.fill_between(
                pos_list, positive_gc_skews, 0, vmin=vmin, vmax=vmax, color="olive"
            )
            gc_skew_track.fill_between(
                pos_list, negative_gc_skews, 0, vmin=vmin, vmax=vmax, color="purple"
            )

            pos_list, telomere_counts = count_special_sequences(seq)
            abs_max_telomere_counts = 50
            vmin, vmax = 0, abs_max_telomere_counts
            telomere_track.fill_between(pos_list, telomere_counts, 0, vmin=vmin, vmax=vmax,color="red")

def concatenate_fasta_records(filepath):
    records = list(SeqIO.parse(filepath, "fasta"))
    concatenated_seq_str = "".join([str(record.seq) for record in records])
    concatenated_seq = Seq(concatenated_seq_str)
    return SeqRecord(concatenated_seq, id=os.path.basename(filepath), description="")

def find_homology_in_folder(folder_path, num_threads):
    print("Starting homology search...")
    fasta_files = [f for f in os.listdir(folder_path) if f.endswith(".fasta")]

    if not os.path.exists("temp"):
        os.mkdir("temp")

    results = []
    lengths = []

    for i in range(len(fasta_files)):
        query = concatenate_fasta_records(os.path.join(folder_path, fasta_files[i]))
        # Store the length and sequence name for the current sequence
        lengths.append({
            "seq": query.id,
            "start": 1,
            "end": len(query.seq)
        })

    def process_pair(pair):
        i, j = pair
        temp_results = []

        # Generating unique file names using thread identifier
        thread_id = threading.get_ident()
        unique_query_filename = f"temp/temp_query_{i}_{thread_id}.fasta"
        unique_subject_filename = f"temp/temp_subject_{j}_{thread_id}.fasta"
        unique_result_filename = f"temp/results_{i}_{j}_{thread_id}.txt"

        query = concatenate_fasta_records(os.path.join(folder_path, fasta_files[i]))
        SeqIO.write(query, unique_query_filename, "fasta")

        subject = concatenate_fasta_records(os.path.join(folder_path, fasta_files[j]))
        SeqIO.write(subject, unique_subject_filename, "fasta")

        blastn_cline = NcbiblastnCommandline(query=unique_query_filename,
                                             subject=unique_subject_filename,
                                             outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                                             out=unique_result_filename, evalue=1e-100)
        blastn_cline()

        with open(unique_result_filename, "r") as result_file:
            for line in result_file:
                parts = line.strip().split("\t")
                q_start, q_end = int(parts[6]), int(parts[7])
                s_start, s_end = int(parts[8]), int(parts[9])
                percent_identity = float(parts[2])
                temp_results.append({
                    "seq1": parts[0],
                    "start1": q_start,
                    "end1": q_end,
                    "seq2": parts[1],
                    "start2": s_start,
                    "end2": s_end,
                    "percent_id": percent_identity
                })

        # Optionally, delete the temporary files after processing
        os.remove(unique_query_filename)
        os.remove(unique_subject_filename)
        os.remove(unique_result_filename)

        return temp_results

    pairs = [(i, j) for i in range(len(fasta_files)) for j in range(i + 1, len(fasta_files))]

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(process_pair, pair) for pair in pairs]
        for future in concurrent.futures.as_completed(futures):
            results.extend(future.result())

    pd.DataFrame(results).to_csv("homology_results.csv", index=False)
    print("Homology search completed.")
    print("Writing sequence lengths...")  # Add this print statement
    pd.DataFrame(lengths).to_csv("sequence_lengths.csv", index=False)
def write_fasta_genome(seqid2seq, outfolder):
    for seq in seqid2seq:
        filename=f"{seq}.fasta"
        seq_record = SeqRecord(seqid2seq[seq], id=seq, description="")
        SeqIO.write(seq_record, os.path.join(outfolder, filename), "fasta")

if __name__ == "__main__":
    args = argparser()
    gbk = Genbank(args.genbank)
    min_size = 10000
    sizes= {key: value for key, value in gbk.get_seqid2size().items() if value >= min_size}
    circos=Circos(sizes, space=1)
    ColorCycler.set_cmap("gist_rainbow")
    chr_names = [s.name for s in circos.sectors]
    colors = ColorCycler.get_color_list(len(chr_names))
    chr_name2color = {name: color for name, color in zip(chr_names, colors)}
    seqid2features = gbk.get_seqid2features(feature_type=None)
    seqid2seq = gbk.get_seqid2seq()

    add_gc_content_skew_tracks(circos, seqid2seq)
    pos, cont=gbk.calc_gc_content()





    for sector in circos.sectors:
        color = chr_name2color[sector.name]
        outer_track = sector.add_track((98, 100))
        outer_track.axis(fc=color)
        major_interval = 1000000
        minor_interval = major_interval/ 10
        if sector.size > major_interval:
            sector.text(sector.name.replace("_", " ").replace("s", "S"), size=10, r=110, color="gray")
            outer_track.xticks_by_interval(major_interval, label_formatter=lambda v: f"{v / 1000000:.0f} Mb")
            outer_track.xticks_by_interval(minor_interval, tick_length=1, show_label=False)
        else:
            outer_track.xticks_by_interval(major_interval, show_label=False)
            outer_track.xticks_by_interval(minor_interval, tick_length=1, show_label=False)
        f_cds_track = sector.add_track((88, 95), r_pad_ratio=0.1)
        r_cds_track = sector.add_track((81, 88), r_pad_ratio=0.1)
        trna_track = sector.add_track((74, 81), r_pad_ratio=0.1)

        for feature in seqid2features[sector.name]:
            if feature.type == "CDS":
                if feature.location.strand == 1:
                    f_cds_track.genomic_features([feature], fc="tomato", )
                else:
                    r_cds_track.genomic_features([feature], fc="skyblue")
            elif feature.type == "tRNA":
                trna_track.genomic_features([feature], color="magenta", lw=0.5)

    tmp_folder = "tmp_circos"
    os.makedirs(tmp_folder, exist_ok=True)
    write_fasta_genome(seqid2seq, tmp_folder)

    find_homology_in_folder(tmp_folder, args.threads)

    with open("homology_results.csv") as f:
        f.readline()
        for line in f:
            line = line.rstrip().split(",")
            name1, start1, end1, name2, start2, end2, percent_id = line
            start1, end1, start2, end2 = map(int, [start1, end1, start2, end2])
            name1, name2 = name1.split(".")[0], name2.split(".")[0]
            coord1=(name1, start1, end1)
            coord2=(name2, start2, end2)
            if len(seqid2seq[name1]) > 100000 and len(seqid2seq[name2]) > 100000 and float(percent_id) > 50 and end1-start1>1000 and end2-start2>1000:
                circos.link(coord1, coord2, color=chr_name2color[name1])



    circos.savefig(args.output if args.output else "circos_plot.svg")
    shutil.rmtree(tmp_folder)
    os.remove("homology_results.csv")
    os.remove("sequence_lengths.csv")
    shutil.rmtree("temp")
