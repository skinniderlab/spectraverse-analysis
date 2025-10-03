import pandas as pd
import numpy as np
import os, sys
from tqdm import tqdm

print("Step2-2: Cleaning metadata and MGF files by removing spectra with identical fragment intensities")

pd.options.mode.chained_assignment = None

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

remove_spectra = []

for i in tqdm(range(len(metadata_index)), desc="Processing metadata"):
    num_peaks = metadata['NUM_PEAKS'][i]
    num_sublists = [sublist for sublist in info[i] if sublist and sublist[0].isdigit()]

    fragment_intens = []
    for sublist in num_sublists:
        fragment_intens.append(float(sublist.split()[1]))

    if all(intens == fragment_intens[0] for intens in fragment_intens) and len(fragment_intens) > 1:
        remove_spectra.append(i)
          
metadata_index_updated = metadata_index.difference(remove_spectra)

info_fin = [info[i] for i in metadata_index_updated]
info_rem = [info[i] for i in remove_spectra]

metadata_cleaned = metadata.drop(remove_spectra).reset_index(drop=True)

with open(output_mgf_dir, 'w') as file:
    for sublist in info_fin:
        for item in sublist:
            file.write("%s\n" % item)
        file.write("\n")    

metadata_cleaned.to_csv(output_csv_dir, index=False)