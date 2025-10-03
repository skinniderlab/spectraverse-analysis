import numpy as np
import os, glob, re
import pandas as pd
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm
from multiprocessing import Pool

print("Step2-7: Removal of spectra with all fragment mass > PRECURSOR_MZ")

input_csv_dir = sys.argv[1]
input_mgf_dir = sys.argv[2]

output_csv_dir = sys.argv[3]
output_mgf_dir = sys.argv[4]

metadata = pd.read_csv(input_csv_dir)

def extract_info(filename):
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
    return info_list

info_before = extract_info(input_mgf_dir)

entries = []
for lst in info_before:
    parent_mass = None
    for item in lst:
        if item.startswith('PRECURSOR_MZ='):
            parent_mass = float(item.split('=')[1])
            break

    peak_mzs = []
    for item in lst:
        parts = item.split()
        if parts and parts[0].replace('.', '', 1).isdigit():
            peak_mzs.append(float(parts[0]))

    if all(mz > parent_mass for mz in peak_mzs):
        entry = {}
        for item in lst:
            if item.startswith('INDEX='):
                entry['INDEX'] = item.split('=', 1)[1]
            if item.startswith('SOURCE='):
                entry['SOURCE'] = item.split('=', 1)[1]
        if 'INDEX' in entry and 'SOURCE' in entry:
            entries.append((entry['INDEX'], entry['SOURCE']))


metadata_select = metadata[~metadata.apply(lambda row: (str(row['INDEX']), str(row['SOURCE'])) in entries, axis=1)]

metadata_index = metadata_select.index

def extract_info_list(filename, indices):
    info_list = []
    record = False
    index = 0
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
                index += 1
            elif record:
                current_info.append(line)
    return [info_list[i] for i in indices]

info = extract_info_list(input_mgf_dir, metadata_index)

metadata_select.to_csv(output_csv_dir, index=False)

with open(output_mgf_dir, 'w') as file:
    for sublist in info:
        for item in sublist:
            file.write("%s\n" % item)
        file.write("\n") 