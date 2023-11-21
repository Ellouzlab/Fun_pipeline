import os
import json
import pandas as pd
from annot_fungi import *

target = "tw_fun_assem"


def read_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data


def get_file(folder_path, ext):
    for file in os.listdir(folder_path):
        if file.endswith(ext):
            return f"{folder_path}/{file}"


def flatten_json(y):
    out = {}

    def flatten(x, name=''):
        if type(x) is dict:
            for a in x:
                flatten(x[a], name + a + '_')
        else:
            out[name[:-1]] = x

    flatten(y)
    return out


# Dictionary to store the data
data = {}

for path in os.listdir(target):
    FUNG = fungi(f"{target}/{path}")
    annotate_result_path = f"{FUNG.predict_dir}/annotate_results"
    json_path = get_file(annotate_result_path, "json")

    # Read and flatten the JSON data
    json_data = read_json(json_path)
    flattened_data = flatten_json(json_data)

    # Iterate over each key in the flattened data and store it
    for key in flattened_data:
        if key not in data:
            data[key] = []
        data[key].append(flattened_data[key])

# Convert the dictionary to a DataFrame
df = pd.DataFrame(data, index=[FUNG.basename for _ in os.listdir(target)])

# Transpose the DataFrame to get basenames as columns
df = df.transpose()

# Print the DataFrame
print(df)

# Print the DataFrame
df.to_csv("annotation_summary.csv")
