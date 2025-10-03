import numpy as np
import os, glob, re
import pandas as pd
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm
from multiprocessing import Pool

print("Step2-4: Removal of unwanted spectra based on various criteria")

input_csv_dir = sys.argv[1]
input_mgf_dir = sys.argv[2]

benchmark_dir = os.path.dirname(input_csv_dir)

temp_mgf_dir = benchmark_dir + '/test-can-temp.mgf'

output_csv_dir = sys.argv[3]
output_mgf_dir = sys.argv[4]

removed_index = []

metadata = pd.read_csv(input_csv_dir)

def convert_to_ms2(value):
    if value == '2' or value == 2:
        return 'MS2'
    return value

metadata['MS_LEVEL'] = metadata['MS_LEVEL'].apply(convert_to_ms2)

index = 0
inside_block = False
smiles_present = False
inchikey_present = False
block_lines = []

with open(input_mgf_dir, 'r') as source, open(temp_mgf_dir, 'w') as destination:
    for line in source:
        if "BEGIN IONS" in line:
            inside_block = True
            block_lines.append(line)
        elif inside_block and "END IONS" in line:
            inside_block = False    
            if not smiles_present:
                block_lines.append(f'SMILES={metadata["CANONICAL_SMILES"][index]}\n')
            if not inchikey_present:
                block_lines.append(f'INCHIKEY={metadata["CANONICAL_INCHIKEY"][index]}\n')     
            block_lines.append(line)
            for block_line in block_lines:
                destination.write(block_line)
            block_lines = []
            index = index + 1
            smiles_present = False
            inchikey_present = False
        elif "SMILES" in line:    
            block_lines.append(f'SMILES={metadata["CANONICAL_SMILES"][index]}\n')   
            smiles_present = True
        elif "INCHIKEY" in line:    
            block_lines.append(f'INCHIKEY={metadata["CANONICAL_INCHIKEY"][index]}\n')   
            inchikey_present = True     
        elif inside_block:
            block_lines.append(line)
        else:
            destination.write(line)


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

metadata = metadata[~metadata['PRECURSOR_MZ'].isnull()]

metadata = metadata[metadata['NUM_PEAKS'] != 0]

metadata = metadata[metadata['MS_LEVEL'] == 'MS2']


remove_neutralloss_adduct = ['[M-H2O]+', '[M-H2O-H]-','[M-2H2O+H]+','[M+H-NH3]+','[M-CH3]-',
                            '[M-2H+H2O]-', '[M-H-C6H10O5}-', '[M-H2O+H]+', '[M-C6H10O5-H]-', '[M+H-H20]', '[M+H-H2O]', 
                            '[M+H-2H2O]', '[M-H2O]', '[M-H20+H]', '[M-H2O-H]', '[M-H2O+H]','[M-2(H2O)+H]', '[M+H-H2O]+', 
                            '[M-H-C6H10O5]', '[M-H-C6H10O5]-', '[M+H-C9H10O5]+', '[M-C6H10O5+H]+', '[M-CO2-H]-', 
                            '[M-H20-H]-', '[M+H-H2O]-', '[M-3H2O+H]+', '[M-4H2O+H]+', '[M-2H2O+NH4]+', '[2M-3H2O+H]+', 
                            '[2M-H2O+H]+', '[M-5H2O+H]+', '[M-2H2O+2H]2+', '[M-3H2O+2H]2+', '[2M-2H2O+H]+']

metadata = metadata[~metadata['ADDUCT'].isin(remove_neutralloss_adduct)]

remove_multimer_adduct = ['[3M-H]-', '[2M+FA-H]-', '[2M-2H2O-H]-', '[2M+H]+', '[2M+Na]+', 'M+2Na', '[2M+Na]', 
                        '[2M+H]', '[2M-H]', '[2M+K]', '[2M-H+Na]', '[M2+Na]', '[M2+H]', '[M+2H+2]',
                        '[M+2]', '[M+2Na]', '[M2+]', '[2M+Ca]2+', 
                        '[2M-H+2Na]+', '[2M+Ca-H]+', '[3M+Ca]2+', '[3M+Na]+', '[3M+Ca-H]+', '[4M+Ca]2+', '[3M+K]+', 
                        '[5M+Ca]2+', '[3M+NH4]+', '[2M+H+CH3CN]+', 
                        '[2M-H]-','[M-H]2-', '[M+2H]2+', '[M+H]2+', '[M+H+Na]2+']

metadata = metadata[~metadata['ADDUCT'].isin(remove_multimer_adduct)]

unknown_adduct = ['[M+H-99]+', '[M+H-99]', '[M-HAc]', '[M+TFA-H]', '[M-CO+H]', '[M+H-(C12H20O9)]+', '[M+Ca]2+', '[M+Ca-H]+', '[M+CH3COO]-/[M-CH3]-']

metadata = metadata[~metadata['ADDUCT'].isin(unknown_adduct)]

def update_ionmode(row):
    adduct = row['ADDUCT']
    ionmode = row['IONMODE']
    if pd.notna(adduct) and adduct.endswith(('+', '-')):
        sign = adduct[-1]
        if sign == '+' and ionmode != 'positive':
            return 'positive'
        elif sign == '-' and ionmode != 'negative':
            return 'negative'
    return ionmode

metadata['IONMODE'] = metadata.apply(update_ionmode, axis=1)

metadata_index = metadata.index

info = extract_info(temp_mgf_dir, metadata_index)

with open(output_mgf_dir, 'w') as file:
    for sublist in info:
        for item in sublist:
            file.write("%s\n" % item)
        file.write("\n") 

metadata.to_csv(output_csv_dir, index=False)