import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
import pandas as pd
import pycircos
import argparse
import concurrent.futures
import threading
from multiprocessing import cpu_count
def arg_parser():
    parser = argparse.ArgumentParser(description="Find homology between sequences in a folder")
    parser.add_argument("-i", "--indir", help="Input directory with fasta files", required=True)
    parser.add_argument("-o", "--output", help="Output directory", required=True)
    parser.add_argument("-t", "--threads", type=int, required=False, default=cpu_count()-1, help="Number of threads to use")
    parser.add_argument('--rerun', action='store_true', help='Rerun homology search without blasting again')
    parser.add_argument('--min_percent_id', type=float, required=False, default=90, help='Minimum percentage identity for homology')
    parser.add_argument('--min_alignment_length', type=int, required=False, default=1000, help='Minimum alignment length for homology')
    return parser.parse_args()
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
        query, query_position_map = concatenate_fasta_records(os.path.join(folder_path, fasta_files[i]))
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

def get_color_gradient(percent_identity, min_percent_identity):
    """
    Returns a color gradient from green (high percent identity) to yellow to red (low percent identity).

    Parameters:
    percent_identity (float): The percent identity value.
    min_percent_identity (float): The minimum percent identity to consider.

    Returns:
    tuple: RGB color.
    """
    # Normalize the percent identity to range between 0 and 1
    normalized_identity = (percent_identity - min_percent_identity) / (100 - min_percent_identity)

    # Invert normalized identity so high values correspond to green and low to red
    inverted_identity = 1 - normalized_identity

    if inverted_identity > 0.5:
        # Closer to red, transition from yellow to red
        red_amount = 1
        green_amount = 2 * (1 - inverted_identity)
    else:
        # Closer to green, transition from green to yellow
        red_amount = 2 * inverted_identity
        green_amount = 1

    return (red_amount, green_amount/1.4, 0)


if __name__ == "__main__":
    args = arg_parser()
    if not args.rerun:
        folder_path = args.indir
        find_homology_in_folder(folder_path, args.threads)

    print('Reading sequence lengths...')
    # Visualization with pycircos
    circle = pycircos.Gcircle(figsize=(9, 10))


    with open("sequence_lengths.csv") as f:
        f.readline()
        for line in f:
            line = line.rstrip().split(",")
            name, length = line[0], int(line[-1])
            arc = pycircos.Garc(arc_id=name, size=length, interspace=2, raxis_range=(935, 985), labelposition=80, label_visible=True)
            circle.add_garc(arc)

    circle.set_garcs(-65, 245)
    for arc_id in circle.garc_dict:
        circle.tickplot(arc_id, raxis_range=(985, 1000), tickinterval=20000, ticklabels=None)

    print('Adding homology results...')
    with open("homology_results.csv") as f:
        f.readline()
        for line in f:
            line = line.rstrip().split(",")
            name1, start1, end1, name2, start2, end2, percent_id = line
            start1, end1, start2, end2 = map(int, [start1, end1, start2, end2])
            percent_id = float(percent_id)
            align_length = max(end1 - start1, end2 - start2)

            if percent_id >= args.min_percent_id and align_length >= args.min_alignment_length:
                source = (name1, start1 - 1, end1, 900)
                destination = (name2, start2 - 1, end2, 900)

                color = get_color_gradient(percent_id, args.min_percent_id)
                circle.chord_plot(source, destination, facecolor=color)

    print('Making Figure...')

    circle.figure.savefig(args.output)
