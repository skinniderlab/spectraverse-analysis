import numpy as np
import os, glob
import pandas as pd
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm
from multiprocessing import Pool

print("Step3-3: Selecting duplicate spectra based on InChIKey first 14 characters and IonMode")

input_csv_dir = sys.argv[1]
input_mgf_dir = sys.argv[2]

mgf_dir = sys.argv[3]
numpy_dir = sys.argv[4]

removed_index = []

metadata = pd.read_csv(input_csv_dir)
metadata_index = metadata.index

def create_dir_and_delete_files(dir_path):
    os.makedirs(dir_path, exist_ok=True)
    files = glob.glob(f'{dir_path}/*')
    for f in files:
        os.remove(f)

create_dir_and_delete_files(mgf_dir)

dedup = metadata.copy()

dedup['INCHI14'] = dedup['INCHIKEY'].str[:14]
dedup = dedup.assign(unique_identifier = dedup[['INCHI14',
                                                'IONMODE']].astype(str).\
                     agg('_'.join, axis=1))

dup = dedup[dedup.duplicated(subset=['unique_identifier'])]
dup_identifiers = dup['unique_identifier'].unique()
dup_df = dedup[dedup['unique_identifier'].isin(dup_identifiers)]

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

def process_identifier(identifier):
    rows = dup_df[dup_df['unique_identifier'] == identifier]
    rows_index = rows.index
    mgf_filepath = mgf_dir + '/{}.mgf'.format(identifier)
    info_temp = [info[i] for i in rows_index]
    assert rows.shape[0] > 1

    with open('{}'.format(mgf_filepath), 'w') as file:
        for sublist in info_temp:
            for item in sublist:
                file.write("%s\n" % item)
            file.write("\n")

with Pool() as pool:
    list(tqdm(pool.imap(process_identifier, dup_identifiers), total=len(dup_identifiers)))
