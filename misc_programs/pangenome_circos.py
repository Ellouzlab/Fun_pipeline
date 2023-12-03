import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
import pandas as pd
import pycircos
import collections
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt


def concatenate_fasta_records(filepath):
    records = list(SeqIO.parse(filepath, "fasta"))
    concatenated_seq_str = "".join([str(record.seq) for record in records])
    concatenated_seq = Seq(concatenated_seq_str)

    position_map = {}
    start = 0
    for record in records:
        position_map[record.id] = (start, start + len(record.seq))
        start += len(record.seq)

    return SeqRecord(concatenated_seq, id=os.path.basename(filepath), description=""), position_map


def find_homology_in_folder(folder_path):
    print("Starting homology search...")
    fasta_files = [f for f in os.listdir(folder_path)]

    if not os.path.exists("temp"):
        os.mkdir("temp")

    results = []
    lengths = []

    for i in range(len(fasta_files)):
        query, query_position_map = concatenate_fasta_records(os.path.join(folder_path, fasta_files[i]))
        SeqIO.write(query, "temp/temp_query.fasta", "fasta")

        # Store the length and sequence name for the current sequence
        lengths.append({
            "seq": query.id,
            "start": 1,
            "end": len(query.seq)
        })

        for j in range(i + 1, len(fasta_files)):
            subject, subject_position_map = concatenate_fasta_records(os.path.join(folder_path, fasta_files[j]))
            SeqIO.write(subject, "temp/temp_subject.fasta", "fasta")

            blastn_cline = NcbiblastnCommandline(query="temp/temp_query.fasta", subject="temp/temp_subject.fasta",
                                                 outfmt=6, out="temp/results.txt", evalue=0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001)
            blastn_cline()

            min_alignment_length = 0  # Set this to your desired minimum alignment length

            with open("temp/results.txt", "r") as result_file:
                for line in result_file:
                    parts = line.strip().split("\t")
                    q_start, q_end = int(parts[6]), int(parts[7])
                    s_start, s_end = int(parts[8]), int(parts[9])
                    alignment_length = max(q_end - q_start + 1, s_end - s_start + 1)  # Calculate alignment length
                    if alignment_length >= min_alignment_length:  # Check if alignment length is above the threshold
                        results.append({
                            "seq1": parts[0],
                            "start1": q_start,
                            "end1": q_end,
                            "seq2": parts[1],
                            "start2": s_start,
                            "end2": s_end
                        })

    # Save the homology results and the sequence lengths to CSV files
    pd.DataFrame(results).to_csv("homology_results.csv", index=False)

    print("Writing sequence lengths...")  # Add this print statement
    pd.DataFrame(lengths).to_csv("sequence_lengths.csv", index=False)


if __name__ == "__main__":
    folder_path = 'pan'
    fasta_files = [f for f in os.listdir(folder_path) if f.endswith(".fasta") or f.endswith(".fa")]
    print(f"Identified {len(fasta_files)} FASTA files: {fasta_files}")  # Add this print statement
    find_homology_in_folder(folder_path)


    Garc = pycircos.Garc
    Gcircle = pycircos.Gcircle
    # Set chromosomes
    circle = Gcircle(figsize=(9, 10))
    with open("sequence_lengths.csv") as f:
        f.readline()
        for line in f:
            line = line.rstrip().split(",")
            name = line[0]
            length = int(line[-1])
            arc = Garc(arc_id=name, size=length, interspace=2, raxis_range=(935, 985), labelposition=80,
                       label_visible=True)
            circle.add_garc(arc)

    circle.set_garcs(-65, 245)
    for arc_id in circle.garc_dict:
        circle.tickplot(arc_id, raxis_range=(985, 1000), tickinterval=20000, ticklabels=None)
    circle.figure

    TRANSPARENCY_LEVEL = 0.5

    values_all = []
    arcdata_dict = collections.defaultdict(dict)
    with open("homology_results.csv") as f:
        f.readline()
        for line in f:
            line = line.rstrip().split(",")
            name1 = line[0]
            start1 = int(line[1]) - 1
            end1 = int(line[2])
            name2 = line[3]
            start2 = int(line[4]) - 1
            end2 = int(line[5])
            source = (name1, start1, end1, 900)
            destination = (name2, start2, end2, 900)

            # Use the facecolor with added transparency for the chord_plot
            facecolor = circle.garc_dict[name1].facecolor
            # If facecolor is a string (e.g., 'blue', '#00FF00'), convert it to RGB format
            if isinstance(facecolor, str):
                facecolor = mcolors.to_rgb(facecolor)

            if isinstance(facecolor, tuple) or isinstance(facecolor, list):
                if len(facecolor) == 3:  # RGB format
                    transparent_color = (*facecolor, TRANSPARENCY_LEVEL)  # Convert to RGBA
                elif len(facecolor) == 4:  # RGBA format
                    transparent_color = (*facecolor[:3], TRANSPARENCY_LEVEL)  # Adjust alpha channel
                else:
                    raise ValueError("Unexpected facecolor format")
            else:
                raise ValueError(f"Unexpected facecolor type: {type(facecolor)}, value: {facecolor}")

            circle.chord_plot(source, destination, facecolor=transparent_color)

    circle.figure
    circle.figure.savefig("tutorial.svg")