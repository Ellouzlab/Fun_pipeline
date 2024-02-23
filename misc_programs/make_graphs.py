import os
import json
import pandas as pd
from annot_fungi import *
from Bio import SeqIO
from pycirclize import Circos
from pycirclize.parser import Genbank
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors


target = "tw_fun_assem"

def find_ext(folder_path, ext):
    for file in os.listdir(folder_path):
        if file.endswith(ext):
            return f"{folder_path}/{file}"


def make_circos_plot(gbk_file, basename):
    gbk = Genbank(gbk_file)

    circos = Circos(sectors={gbk.name: gbk.range_size})
    circos.text(basename, size=12, r=20)
    sector = circos.get_sector(gbk.name)

    major_ticks_interval = 10000000
    minor_ticks_interval = 500000
    outer_track = sector.add_track((98, 100))
    outer_track.axis(fc="lightgrey")
    outer_track.xticks_by_interval(
        major_ticks_interval, label_formatter=lambda v: f"{v / 10 ** 6:.1f} Mb"  # Adjusted label size
    )
    outer_track.xticks_by_interval(minor_ticks_interval, tick_length=1, show_label=False)

    # Plot Forward CDS, Reverse CDS, rRNA, tRNA
    f_cds_track = sector.add_track((90, 97), r_pad_ratio=0.1)
    f_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=1), fc="red")

    r_cds_track = sector.add_track((83, 90), r_pad_ratio=0.1)
    r_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=-1), fc="blue")

    trna_track = sector.add_track((76, 83), r_pad_ratio=0.1)
    trna_track.genomic_features(gbk.extract_features("tRNA"), color="magenta", lw=0.1)

    # Plot GC content
    gc_content_track = sector.add_track((50, 65))

    pos_list, gc_contents = gbk.calc_gc_content()
    gc_contents = gc_contents - gbk.calc_genome_gc_content()
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

    # Plot GC skew
    gc_skew_track = sector.add_track((35, 50))

    pos_list, gc_skews = gbk.calc_gc_skew()
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
    outpath = f"Results/circos/{basename}.png"

    fig = circos.plotfig()
    handles = [
        Patch(color="red", label="Forward CDS"),
        Patch(color="blue", label="Reverse CDS"),
        Patch(color="green", label="rRNA"),
        Patch(color="magenta", label="tRNA"),
        Line2D([], [], color="black", label="Positive GC Content", marker="^", ms=6, ls="None"),
        Line2D([], [], color="grey", label="Negative GC Content", marker="v", ms=6, ls="None"),
        Line2D([], [], color="olive", label="Positive GC Skew", marker="^", ms=6, ls="None"),
        Line2D([], [], color="purple", label="Negative GC Skew", marker="v", ms=6, ls="None"),
    ]
    _ = circos.ax.legend(handles=handles, bbox_to_anchor=(0.5, 0.475), loc="center", fontsize=8)
    circos.savefig(outpath, dpi=300)

    # Add legend

    print(f"{gbk_file}.png")

def parse_gbk_for_cog(gbk_file, basename):
    cog_counts = {}

    with open(gbk_file, 'r') as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    note = feature.qualifiers.get('note', [])
                    for item in note:
                        if "COG:" in item:
                            cogs = item.split("COG:")[1].split(';')[0]
                            for cog in cogs.split(','):
                                cog = cog.strip()
                                cog_counts[cog] = cog_counts.get(cog, 0) + 1

    # Create a DataFrame with basename as the index and COG categories as columns
    df = pd.DataFrame([cog_counts], index=[basename])

    return df


def get_unique_colors(num_colors):
    """Generate a list of distinct colors."""
    colors = plt.cm.get_cmap('tab20', num_colors)  # Use 'tab20' as the base, which has 20 unique colors
    return colors(np.linspace(0, 1, num_colors))

def get_extended_unique_colors(num_colors):
    """Generate a list of distinct colors from multiple colormaps."""
    base_colormaps = ['tab20', 'tab20b', 'tab20c', 'Set3', 'Pastel1']  # List of colormaps
    colors = []

    for cmap_name in base_colormaps:
        cmap = plt.cm.get_cmap(cmap_name)
        colors.extend(cmap(np.linspace(0, 1, cmap.N)))

    if num_colors > len(colors):
        # If more colors are needed, add random colors
        np.random.seed(0)  # For reproducibility
        additional_colors = np.random.rand(num_colors - len(colors), 3)
        colors.extend(additional_colors)

    return colors[:num_colors]

def plot_stacked_bar_percentage(df):
    # Convert gene counts to percentages

    cog_descriptions = {
        "A": "A: RNA processing and modification",
        "B": "B: Chromatin structure and dynamics",
        "C": "C: Energy production and conversion",
        "D": "D: Cell cycle control and mitosis",
        "E": "E: Amino acid transport and metabolism",
        "F": "F: Nucleotide transport and metabolism",
        "G": "G: Carbohydrate transport and metabolism",
        "H": "H: Coenzyme transport and metabolism",
        "I": "I: Lipid transport and metabolism",
        "J": "J: Translation, ribosomal structure and biogenesis",
        "K": "K: Transcription",
        "L": "L: Replication and repair",
        "M": "M: Cell wall/membrane biogenesis",
        "N": "N: Cell motility",
        "O": "O: Post-translational modification, protein turnover, chaperones",
        "P": "P: Inorganic ion transport and metabolism",
        "Q": "Q: Secondary metabolites biosynthesis, transport, and catabolism",
        "R": "R: General function prediction only",
        "S": "S: Function unknown",
        "T": "T: Signal transduction mechanisms",
        "U": "U: Intracellular trafficking and secretion",
        "V": "V: Defense mechanisms",
        "W": "W: Extracellular structures",
        "Y": "Y: Nuclear structure",
        "Z": "Z: Cytoskeleton",
        # Add other categories if you have more
    }
    # Convert gene counts to percentages
    # Convert gene counts to percentages
    df_percentage = df.div(df.sum(axis=1), axis=0) * 100

    # Get unique colors for the number of categories
    num_categories = len(df.columns)
    unique_colors = get_extended_unique_colors(num_categories)

    # Create a wider stacked bar plot with custom colors
    ax = df_percentage.plot(kind='bar', stacked=True, figsize=(20, 8), width=0.8, color=unique_colors)

    # Set up the line graph for total counts
    ax2 = ax.twinx()
    total_count = df.sum(axis=1)
    total_count.plot(kind='line', marker='o', color='black', ax=ax2, linewidth=2)
    ax2.set_ylabel("Total COG Annotations", fontsize=16)
    ax2.tick_params(axis='y', labelsize=14)

    # Set both y-axes to start from zero
    ax.set_ylim(bottom=0)
    ax2.set_ylim(bottom=0)

    # Set font size for all plot elements
    plt.rcParams.update({'font.size': 14})

    # Adding labels and title with larger font size
    ax.set_xlabel("Fungal Genomes", fontsize=16)
    ax.set_ylabel("Percentage of Genes (%)", fontsize=16)

    # Format y-axis to show percentage
    ax.yaxis.set_major_formatter(mticker.PercentFormatter())

    # Combine handles and labels for both bar and line plots for the legend
    handles, labels = ax.get_legend_handles_labels()
    line_label = "Total COG Annotations"
    line_handle, = ax2.get_lines()
    handles.append(line_handle)
    labels.append(line_label)

    # Sort legend entries to match the order of the bars and include line plot
    ordered_labels_handles = sorted(zip(labels, handles), key=lambda lh: df_percentage.columns.tolist().index(lh[0]) if lh[0] in df_percentage.columns else -1, reverse=True)
    ordered_handles = [handle for label, handle in ordered_labels_handles]
    ordered_labels = [cog_descriptions.get(label, label) for label, handle in ordered_labels_handles]

    # Create a smaller legend with smaller font size and adjust position
    ax.legend(ordered_handles, ordered_labels, bbox_to_anchor=(1.2, 1), loc='upper left', fontsize=12, title_fontsize=14)

    # Beautify the x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    # Thicken the axes lines
    [i.set_linewidth(1.5) for i in ax.spines.values()]

    # Adjust layout to prevent overlap and ensure everything fits
    plt.tight_layout()

    # Display the plot
    plt.savefig("test_cog.png", dpi=300)

def plot_stacked_bar_percentage_go(df, type):
    # Convert gene counts to percentages
    df_percentage = df.div(df.sum(axis=1), axis=0) * 100

    # Get unique colors for the number of categories
    num_categories = len(df.columns)
    unique_colors = get_extended_unique_colors(num_categories)

    # Create a wider stacked bar plot with custom colors
    ax = df_percentage.plot(kind='bar', stacked=True, figsize=(20, 8), width=0.8, color=unique_colors)

    # Set up the line graph for total counts
    ax2 = ax.twinx()
    total_count = df.sum(axis=1)
    total_count.plot(kind='line', marker='o', color='black', ax=ax2, linewidth=2)
    ax2.set_ylabel("Total Go Process Annotations", fontsize=16)
    ax2.tick_params(axis='y', labelsize=14)

    # Set both y-axes to start from zero
    ax.set_ylim(bottom=0)
    ax2.set_ylim(bottom=0)

    # Set font size for all plot elements
    plt.rcParams.update({'font.size': 14})

    # Adding labels and title with larger font size
    ax.set_xlabel("Samples", fontsize=16)
    ax.set_ylabel("Percentage of Genes (%)", fontsize=16)

    # Format y-axis to show percentage
    ax.yaxis.set_major_formatter(mticker.PercentFormatter())

    # Combine handles and labels for both bar and line plots for the legend
    handles, labels = ax.get_legend_handles_labels()
    line_label = "Total Annotations"
    line_handle, = ax2.get_lines()
    handles.append(line_handle)
    labels.append(line_label)

    # Sort legend entries to match the order of the bars and include line plot
    ordered_labels_handles = sorted(zip(labels, handles), key=lambda lh: df_percentage.columns.tolist().index(lh[0]) if lh[0] in df_percentage.columns else -1, reverse=True)
    ordered_handles = [handle for label, handle in ordered_labels_handles]
    ordered_labels = [label for label, handle in ordered_labels_handles]

    # Create a smaller legend with smaller font size and adjust position
    ax.legend(ordered_handles, ordered_labels, bbox_to_anchor=(1.2, 1), loc='upper left', fontsize=12, title_fontsize=14)

    # Beautify the x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    # Thicken the axes lines
    [i.set_linewidth(1.5) for i in ax.spines.values()]

    # Adjust layout to prevent overlap and ensure everything fits
    plt.tight_layout()

    # Display the plot
    plt.savefig(f"test_annotations_go_{type}.svg", dpi=300)
def parse_gbk_for_GO_process(gbk_file, basename):
    go_counts = {}

    with open(gbk_file, 'r') as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    note = feature.qualifiers.get('note', [])
                    for item in note:
                        item_list=item.split(';')
                        for entry in item_list:
                            if "GO_process" in entry:
                                stuff=entry.split('- ')[1].split(' [')[0]
                                if not stuff in go_counts:
                                    go_counts[stuff] = 1
                                else:
                                    go_counts[stuff] += 1
    sorted_dict = dict(sorted(go_counts.items(), key=lambda item: item[1], reverse=True))

    # Split the sorted dictionary into top 10 and the rest
    top_10_items = dict(list(sorted_dict.items())[:10])
    other_items_sum = sum(list(sorted_dict.values())[10:])

    # Create a new dictionary with top 10 items and an 'other' category
    first_10_elements_with_other = top_10_items
    if other_items_sum > 0:
        first_10_elements_with_other['Other Go Processes'] = other_items_sum

    # Create a DataFrame with basename as the index and GO categories as columns
    df = pd.DataFrame([first_10_elements_with_other], index=[basename])

    return df


def parse_gbk_for_GO(gbk_file, basename, type):
    go_counts = {}

    with open(gbk_file, 'r') as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    note = feature.qualifiers.get('note', [])
                    for item in note:
                        item_list=item.split(';')
                        for entry in item_list:
                            if f"GO_{type}" in entry:
                                stuff=entry.split('- ')[1].split(' [')[0]
                                if not stuff in go_counts:
                                    go_counts[stuff] = 1
                                else:
                                    go_counts[stuff] += 1
    sorted_dict = dict(sorted(go_counts.items(), key=lambda item: item[1], reverse=True))

    # Split the sorted dictionary into top 10 and the rest
    top_10_items = dict(list(sorted_dict.items())[:])


    # Create a DataFrame with basename as the index and GO categories as columns
    df = pd.DataFrame([top_10_items], index=[basename])

    return df


def top_x_columns_with_others(df, x):
    if x >= df.shape[1]:
        print(
            "Warning: x is larger than or equal to the number of columns in the DataFrame. Returning the original DataFrame.")
        return df

    column_sums = df.sum()
    sorted_columns = column_sums.sort_values(ascending=False).index.tolist()
    top_columns = sorted_columns[:x]
    other_columns = sorted_columns[x:]
    df['Other'] = df[other_columns].sum(axis=1)
    return df[top_columns + ['Other']]


def main():
    row_list=[]
    row_list2=[]
    row_list3=[]
    row_list4=[]
    for path in os.listdir(target):
        FUNG = fungi(f"{target}/{path}")

        result_dir_path = f"{FUNG.predict_dir}/annotate_results"

        gbk_path= find_ext(result_dir_path, "gbk")
        name=''.join(('\n'.join(FUNG.basename.split('_'))).split('strain\n'))

        basename=' '.join(FUNG.basename.split("_"))

        #make_circos_plot(gbk_path, name)

        #df = parse_gbk_for_cog(gbk_path, basename)
        #row_list.append(df)

        df2 = parse_gbk_for_GO(gbk_path, basename, 'process')
        row_list2.append(df2)

        df3 = parse_gbk_for_GO(gbk_path, basename, 'component')
        row_list3.append(df3)

        df4 = parse_gbk_for_GO(gbk_path, basename, 'function')
        row_list4.append(df4)

    #df = pd.concat(row_list)
    df2 = pd.concat(row_list2)
    df2 = top_x_columns_with_others(df2, 15)
    df2.to_csv('test_go.csv')
    df3 = pd.concat(row_list3)
    df3 = top_x_columns_with_others(df3, 15)
    df3.to_csv('test_go_component.csv')
    df4 = pd.concat(row_list4)
    df4 = top_x_columns_with_others(df4, 15)
    df4.to_csv('test_go_function.csv')

    #plot_stacked_bar_percentage(df)
    plot_stacked_bar_percentage_go(df2, 'process')
    plot_stacked_bar_percentage_go(df3, 'component')
    plot_stacked_bar_percentage_go(df4, 'function')

if __name__ == "__main__":
    main()