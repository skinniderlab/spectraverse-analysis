import numpy as np
import os, glob, re
import pandas as pd
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm

print("Step2-6: Removal of spectra based on precursor and fragment mass check")

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

for index, row in metadata.iterrows():
    if int(row['NUM_PEAKS']) == 1:
        precursor_mz = row['PRECURSOR_MZ']
        filtered_sublists = [sublist for sublist in info[index] if sublist and sublist[0].isdigit()]
        assert len(filtered_sublists) == 1
        fragment_mass = float(filtered_sublists[0].split()[0])
        if fragment_mass <= 10:
            raise ValueError('Spectra with fragment with mass less than 10 found:', row)
        elif fragment_mass >= 1000:
            raise ValueError('Spectra with fragment with mass greater than 1000 found:', row)   
        min_mz = precursor_mz - 1.6
        if fragment_mass >= min_mz:   
            removed_index.append(index)   
    else:
        filtered_sublists = [sublist for sublist in info[index] if sublist and sublist[0].isdigit()]
        frag_mass_list = []
        for sublist in filtered_sublists:
            frag_mass_list.append(float(sublist.split()[0]))
        if any(frag_mass >= 1000 for frag_mass in frag_mass_list):
            raise ValueError('Spectra with fragment with mass greater than 1000 found:', frag_mass_list)

metadata = metadata[~metadata.index.isin(removed_index)]
metadata_index = metadata.index
info = extract_info(input_mgf_dir, metadata_index)

with open(output_mgf_dir, 'w') as file:
    for sublist in info:
        for item in sublist:
            file.write("%s\n" % item)
        file.write("\n") 

metadata.to_csv(output_csv_dir, index=False)      