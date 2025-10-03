import pandas as pd
import numpy as np
import re, os, sys, math
from rdkit import Chem

print("Step2-8: CSV to MGF conversion (accomplished with metadata modification)")

csv_dir = sys.argv[1]
in_mgf_dir = sys.argv[2]
out_mgf_dir = sys.argv[3]

metadata = pd.read_csv(csv_dir)

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

info = extract_info(in_mgf_dir)

with open(out_mgf_dir, 'w') as mgf_file:
    for i in range(len(metadata)):
        mgf_file.write("BEGIN IONS\n")

        mgf_file.write("INDEX={}\n".format(metadata['INDEX'][i]))
        mgf_file.write("SOURCE={}\n".format(metadata['SOURCE'][i]))
        mgf_file.write("SMILES={}\n".format(metadata['CANONICAL_SMILES'][i].split()[0]))

        ionmode = metadata['IONMODE'][i]
        if not pd.isna(ionmode):
            mgf_file.write("IONMODE={}\n".format(ionmode))

        adduct = metadata['ADDUCT'][i]
        if not pd.isna(adduct):
            mgf_file.write("ADDUCT={}\n".format(adduct))

        # compound_name = metadata['COMPOUND_NAME'][i]
        # if not pd.isna(compound_name):
        #     mgf_file.write("COMPOUND_NAME={}\n".format(compound_name))     
            
        mgf_file.write("NUM_PEAKS={}\n".format(metadata['NUM_PEAKS'][i]))
        mgf_file.write("PRECURSOR_MZ={}\n".format(metadata['PRECURSOR_MZ'][i]))
        mgf_file.write("PARENT_MASS={}\n".format(metadata['EXACT_MASS_CANONICAL'][i]))
        mgf_file.write("MS_LEVEL={}\n".format(metadata['MS_LEVEL'][i]))

        for line in info[i]:
            if re.match(r'^\d', line):
                mgf_file.write(line+'\n')

        mgf_file.write("END IONS\n\n")        
