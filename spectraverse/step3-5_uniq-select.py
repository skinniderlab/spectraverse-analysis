import numpy as np
import os, glob
import pandas as pd
import spectral_entropy
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm
from multiprocessing import Pool

print("Step3-5: Selecting unique spectra from duplicate spectra based on pairwise cosine scores")

input_csv_dir = sys.argv[1]
input_mgf_dir = sys.argv[2]

mgf_dir = sys.argv[3]
numpy_dir = sys.argv[4]

benchmark_dir = os.path.dirname(input_csv_dir)

temp_mgf_dir = benchmark_dir + '/metadata-can-rem-lowres-prec-highfragmass-matchms-tempdot.mgf'

output_csv_dir = sys.argv[5]
output_mgf_dir = sys.argv[6]

removed_index = []

metadata = pd.read_csv(input_csv_dir)
metadata_index = metadata.index

dedup = metadata.copy()

dedup['INCHI14'] = dedup['INCHIKEY'].str[:14]
dedup = dedup.assign(unique_identifier = dedup[['INCHI14',
                                                'IONMODE']].astype(str).\
                     agg('_'.join, axis=1))

dup = dedup[dedup.duplicated(subset=['unique_identifier'])]
dup_identifiers = dup['unique_identifier'].unique()
dup_df = dedup[dedup['unique_identifier'].isin(dup_identifiers)]
dup_df['NUM_PEAKS'] = pd.to_numeric(dup_df['NUM_PEAKS'])

# remove duplicate spectra above a minimum dot-product, keeping the one with the most peaks
min_dot_prod = 0.99

commecial_data = ['nist_23_hr_msms_new.mgf', 'nist_23_hr_msms2_new.mgf', 'agilent_fold_csv_matchms.mgf']
 
def process_identifier(identifier):
    dup_entries = dup_df[dup_df['unique_identifier'] == identifier]
    numpy_filepath = numpy_dir + '/{}.npy'.format(identifier)
    dot_products = np.load(numpy_filepath)
    dot_products_select =  np.triu(dot_products, k=1)
    indices = np.where(dot_products_select > min_dot_prod)
    index_pairs = [list(index) for index in zip(*indices)]
    
    drop_indexes = []
    for pairs in index_pairs:
        row = dup_entries.iloc[pairs]
        if row['SOURCE'].isin(commecial_data).any():
            drop_index = row[row['SOURCE'].isin(commecial_data)].index
            if len(drop_index) > 1:
                drop_index = row['NUM_PEAKS'].idxmin()
        drop_index = row['NUM_PEAKS'].idxmin()
        drop_indexes.append(drop_index)
    drop_indexes = list(set(drop_indexes))
    return drop_indexes

with Pool() as pool:
    results = list(tqdm(pool.imap(process_identifier, dup_identifiers), total=len(dup_identifiers)))

drop_indexes = [index for sublist in results for index in sublist]

dedup.drop(drop_indexes, inplace=True)         


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


dedup_index = dedup.index
info = extract_info(input_mgf_dir, dedup_index)


with open(temp_mgf_dir, 'w') as file:
    for sublist in info:
        for item in sublist:
            file.write("%s\n" % item)
        file.write("\n") 


dedup.to_csv(output_csv_dir, index=False)  

metadata = pd.read_csv(output_csv_dir)

def modify_adduct_lines(input_filename, output_filename, metadata, indices):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        index = 0
        for line in infile:
            if 'ADDUCT=' in line:
                line = f'ADDUCT={metadata.iloc[indices[index]]["ADDUCT"]}\n'
            if 'PARENT_MASS=' in line:
                line = f'PARENT_MASS={metadata.iloc[indices[index]]["PARENT_MASS"]}\n'  
            if 'SMILES=' in line:
                line = f'SMILES={metadata.iloc[indices[index]]["SMILES"]}\n'  
            if 'INCHIKEY=' in line:
                line = f'INCHIKEY={metadata.iloc[indices[index]]["INCHIKEY"]}\n'    
            if 'INCHI=' in line:
                line = f'INCHI={metadata.iloc[indices[index]]["INCHI"]}\n'
            if 'FORMULA=' in line:
                line = f'FORMULA={metadata.iloc[indices[index]]["FORMULA"]}\n'
            if 'CHARGE=' in line:
                line = f'CHARGE={metadata.iloc[indices[index]]["CHARGE"]}\n'            
            if 'END IONS' in line:  
                index += 1
            outfile.write(line)

metadata_index = metadata.index
modify_adduct_lines(temp_mgf_dir, output_mgf_dir, metadata, metadata_index)
