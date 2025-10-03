import numpy as np
import os, glob, re
import pandas as pd
import spectral_entropy
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm
from multiprocessing import Pool

print("Step3-6: Noise removal and checking for spectra with single fragment mass too close to precursor m/z")

input_csv_dir = sys.argv[1]
input_mgf_dir = sys.argv[2]

benchmark_dir = os.path.dirname(input_csv_dir)

temp_mgf_dir = benchmark_dir + '/metadata-can-rem-lowres-prec-highfragmass-matchms-dotcheck-temp.mgf'

output_csv_dir = sys.argv[3]
output_mgf_dir = sys.argv[4]


metadata = pd.read_csv(input_csv_dir)
metadata_index = metadata.index
metadata['NUM_PEAKS_OLD'] = metadata['NUM_PEAKS']

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

info = extract_info(input_mgf_dir)

noise_threshold = 0.001

noiserem_info = []
for i in range(len(info)):
    num_sublists = [sublist for sublist in info[i] if sublist and sublist[0].isdigit()]
    fragment_masses = []
    for sublist in num_sublists:
        fragment_masses.append(float(sublist.split()[0]))
    delta = np.diff(fragment_masses)
    median_delta = np.median(delta)

    num_sublists_noiserem = [
        sublist for sublist in num_sublists
        if not sublist[0].isdigit() or (sublist[0].isdigit() and float(sublist.split()[1]) >= noise_threshold)
    ]

    mz_intensity = np.array([
        [float(sublist.split()[0]), float(sublist.split()[1])]
        for sublist in num_sublists_noiserem
    ])

    if mz_intensity.shape[0] > 4096:
        mz_intensity = mz_intensity[mz_intensity[:, 1].argsort()[::-1]]
        mz_intensity = mz_intensity[:4096]
        mz_intensity = mz_intensity[mz_intensity[:, 0].argsort()]

    num_sublists_noiserem = [
        f"{row[0]} {row[1]}" for row in mz_intensity
    ]

    updated_num_peaks = len(num_sublists_noiserem)

    metadata.loc[i, 'NUM_PEAKS'] = updated_num_peaks

    main_sublists = [sublist for sublist in info[i] if sublist and not sublist[0].isdigit()]

    def update_num_peaks(spectra_list, new_value):
        updated_list = []
        for item in spectra_list:
            if item.startswith('NUM_PEAKS='):
                old_value = item.split('=')[1]
                updated_list.append(f'NUM_PEAKS_OLD={old_value}') 
                updated_list.append(f'NUM_PEAKS={new_value}')
            else:
                updated_list.append(item)
        return updated_list

    updated_main_sublists = update_num_peaks(main_sublists, updated_num_peaks) 

    noiserem_sublists = [sublist for sublist in updated_main_sublists if sublist and not sublist[0].isdigit()]
    noiserem_spectra_list = noiserem_sublists[:-1] + num_sublists_noiserem + noiserem_sublists[-1:]

    noiserem_info.append(noiserem_spectra_list)

assert len(noiserem_info) == len(metadata), "Length of noiserem_info does not match length of metadata."

with open(temp_mgf_dir, 'w') as file:
    for sublist in noiserem_info:
        for item in sublist:
            file.write("%s\n" % item)
        file.write("\n")


info_temp = extract_info(temp_mgf_dir)

removed_index = []

for index, row in metadata.iterrows():
    if int(row['NUM_PEAKS']) == 1:
        precursor_mz = row['PRECURSOR_MZ']
        filtered_sublists = [sublist for sublist in info_temp[index] if sublist and sublist[0].isdigit()]
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
        filtered_sublists = [sublist for sublist in info_temp[index] if sublist and sublist[0].isdigit()]
        frag_mass_list = []
        for sublist in filtered_sublists:
            frag_mass_list.append(float(sublist.split()[0]))
        if any(frag_mass >= 1000 for frag_mass in frag_mass_list):
            raise ValueError('Spectra with fragment with mass greater than 1000 found:', frag_mass_list)

info_temp_removed = [info_temp[i] for i in removed_index]

def extract_info_indexes(filename, indices):
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


metadata = metadata[~metadata.index.isin(removed_index)]
metadata_index = metadata.index
info_out = extract_info_indexes(temp_mgf_dir, metadata_index)

assert len(info_out) == len(metadata), "Length of info_out does not match length of metadata."

metadata.to_csv(output_csv_dir, index=False)

with open(output_mgf_dir, 'w') as file:
    for sublist in info_out:
        for item in sublist:
            file.write("%s\n" % item)
        file.write("\n")