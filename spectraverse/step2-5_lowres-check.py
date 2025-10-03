import numpy as np
import os, glob, re
import pandas as pd
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm

print("Step2-5: Removal of low-resolution spectra")

input_csv_dir = sys.argv[1]
input_mgf_dir = sys.argv[2]

output_csv_dir = sys.argv[3]
output_mgf_dir = sys.argv[4]

metadata = pd.read_csv(input_csv_dir)
metadata_index = metadata.index

def extract_info(filename, indices):
    info_list = []
    record = False
    with open(filename, 'r') as file:
        current_info = []
        for line in file:
            line = line.strip()
            if 'BEGIN IONS' in line:
                record = True
                current_info = [line]
            elif 'END IONS' in line:
                current_info.append(line)
                record = False
                info_list.append(current_info)
            elif record:
                current_info.append(line)
    return [info_list[i] for i in indices]

info = extract_info(input_mgf_dir, metadata_index)

removed_index = []

def get_precision(value):
    value_str = str(value)
    if '.' in value_str:
        decimal_part = value_str.split('.')[1]
        return len(decimal_part)
    else:
        return 0

for index, row in metadata.iterrows():
    prec_lowres = False
    frag_lowres = False

    precursor_mz = row['PRECURSOR_MZ']
    precursor_precision = get_precision(precursor_mz)

    if precursor_precision < 2:
        prec_lowres = True

    filtered_sublists = [sublist for sublist in info[index] if sublist and sublist[0].isdigit()]

    frag_res_list = []
    for sublist in filtered_sublists:
        fragment_mass = float(sublist.split()[0])
        fragment_precision = get_precision(fragment_mass)
        if fragment_precision < 2:
            frag_res_list.append(True)
        else:
            frag_res_list.append(False)    
    if all(frag_res_list):
        frag_lowres = True

    if frag_lowres:
        removed_index.append(index)  


info_removed = [info[i] for i in removed_index]

metadata = metadata[~metadata.index.isin(removed_index)]
metadata_index = metadata.index
info = extract_info(input_mgf_dir, metadata_index)

float_pattern = re.compile(r"[-+]?\d*\.\d+|\d+")

with open(output_mgf_dir, 'w') as file:
    for sublist in info:
        for item in sublist:
            file.write("%s\n" % item)
        file.write("\n") 

metadata.to_csv(output_csv_dir, index=False)