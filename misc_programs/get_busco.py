import os
import json
import pandas as pd
from annot_fungi import *

def grab_string(string_list, key):
    for string in string_list:
        if key in string:
            return string.split(' ')[0]


def read_txt(file):
    with open(file, "r") as f:
        return f.read()

target = "tw_fun_assem"

if __name__ == "__main__":
    row_list = []
    for path in os.listdir(target):
        FUNG = fungi(f"{target}/{path}")

        busco_log=f"{FUNG.log_dir}/busco.log"
        chop=read_txt(busco_log).split('Results:')[1]
        short_chop=chop.split('INFO	')

        complete_BUSCO=int(grab_string(short_chop, 'Complete BUSCOs'))
        complete_single=int(grab_string(short_chop, '(S)'))
        complete_multiple=int(grab_string(short_chop, '(D)'))
        complete_fragmented=int(grab_string(short_chop, '(F)'))
        missing_BUSCO=int(grab_string(short_chop, '(M)'))

        result = pd.DataFrame({'Complete BUSCOs': [complete_BUSCO], 'Complete single-copy BUSCOs': [complete_single], 'Complete duplicated BUSCOs': [complete_multiple], 'Fragmented BUSCOs': [complete_fragmented], 'Missing BUSCOs': [missing_BUSCO]}, index=[' '.join(FUNG.basename.split('_'))])
        row_list.append(result)

    merged_df = pd.concat(row_list)
    merged_df.to_csv("busco.csv")

