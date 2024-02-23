import pandas as pd
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Plot stacked bar graph of metabolite counts")
parser.add_argument("-i", "--input", help="Input file", required=True)
args = parser.parse_args()


df = pd.read_csv(args.input, sep='\t')

# Set the index to the first column to use as labels for the x-axis
df.set_index(df.columns[0], inplace=True)

# Define a new custom color palette that is visually appealing
color_list = [
    '#1f77b4',  # Muted blue
    '#ff7f0e',  # Safety orange
    '#2ca02c',  # Cooked asparagus green
    '#d62728',  # Brick red
    '#9467bd',  # Muted purple
    '#8c564b',  # Chestnut brown
    '#e377c2',  # Raspberry yogurt pink
    '#7f7f7f',  # Middle gray
    '#bcbd22',  # Curry yellow-green
    '#17becf'   # Blue-teal
]

# Assigning colors to each column
color_dict = {df.columns[i]: color_list[i] for i in range(len(df.columns))}

# Plotting the stacked bar chart with custom colors
ax = df.plot(kind='bar', stacked=True, figsize=(14, 8), color=[color_dict.get(x) for x in df.columns])

# Set the labels and title
ax.set_xlabel('Strains', fontsize=12)
ax.set_ylabel('Counts', fontsize=12)
ax.set_title('Stacked Bar Graph of Different Metabolites', fontsize=14)

# Rotate the x-axis labels for better visibility
plt.xticks(rotation=45, ha='right')

# Adjust the legend to match the order of the bars
handles, labels = ax.get_legend_handles_labels()
# Note: The order of handles/labels is reversed to match the order of the stacked bars
ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1, 1), title='Metabolites')

# Tight layout to adjust for the legend
plt.tight_layout()

# Save the plot to a file
plt.savefig('antismash_custom_colors.png')