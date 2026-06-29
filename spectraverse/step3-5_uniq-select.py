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
 
library_order = [
        'tony_AurelieRoux2012_neg_processed.mgf','tony_AurelieRoux2012_pos_processed.mgf',
        'tony_DamienJimenez2019_neg.mgf', 'tony_DamienJimenez2019_pos.mgf',
        'tony_GUIA_neg_172_spectra_processed.mgf','tony_GUIA_pos_449_spectra_processed.mgf',
        'tony_MeRgeIon_processed.mgf','tony_mz_vault_process.mgf',
        'tony_Narayanaswamy2020_MetaboKit_processed.mgf',
        'tony_PFAS-identified_processed.mgf','tony_PFAS-library_processed.mgf',
        'tony_PrasadPhapale2021_Curatr_neg_processed.mgf','tony_PrasadPhapale2021_Curatr_pos_processed.mgf',
        'tony_RomanPopov2023-processed.mgf', 'tony_SMID-DB-processed.mgf',
        'tony_StephanBeisken2014_neg_processed.mgf','tony_StephanBeisken2014_pos_processed.mgf',
        'tony_Zheng2024_neg_processed.mgf','tony_Zheng2024_pos_processed.mgf',

        'tony_FOODB_processed.mgf','tony_HMDB_processed.mgf','tony_MIME_processed.mgf', 

        'tony_20231031_nihnp_library_neg_all_lib_MS2_processed.mgf',
        'tony_20231031_nihnp_library_pos_all_lib_MS2_processed.mgf',
        'tony_20231130_mcescaf_library_neg_all_lib_MS2_processed.mgf',
        'tony_20231130_mcescaf_library_pos_all_lib_MS2_processed.mgf',
        'tony_20231130_otavapep_library_neg_all_lib_MS2_processed.mgf',
        'tony_20231130_otavapep_library_pos_all_lib_MS2_processed.mgf',
        'tony_20240411_mcebio_library_neg_all_lib_MS2_processed.mgf',
        'tony_20240411_mcebio_library_pos_all_lib_MS2_processed.mgf',

        'tony_BioMSMS-Pos-PlaSMA-processed.mgf',
        'tony_BioMSMS-neg-PlaSMA-processed.mgf',
        'tony_KI-GIAR_zic-HILIC_Pos_v0.90-processed.mgf',
        'tony_MSMS-Neg-CASMI2016-processed.mgf',
        'tony_MSMS-Neg-FiehnHILIC-processed.mgf',
        'tony_MSMS-Neg-GNPS-processed.mgf',
        'tony_MSMS-Neg-MassBank-processed.mgf',
        'tony_MSMS-Neg-MassBankEU-processed.mgf',
        'tony_MSMS-Neg-MetaboBASE-processed.mgf',
        'tony_MSMS-Neg-PlaSMA-processed.mgf',
        'tony_MSMS-Neg-Respect-processed.mgf',
        'tony_MSMS-Neg-RikenOxPLs-processed.mgf',
        'tony_MSMS-Neg-Vaniya-Fiehn_Natural_Products_Library_20200109-processed.mgf',
        'tony_MSMS-Pos-CASMI2016-processed.mgf',
        'tony_MSMS-Pos-FiehnHILIC-processed.mgf',
        'tony_MSMS-Pos-GNPS-processed.mgf',
        'tony_MSMS-Pos-MassBank-processed.mgf',
        'tony_MSMS-Pos-MassBankEU-processed.mgf',
        'tony_MSMS-Pos-MetaboBASE-processed.mgf',
        'tony_MSMS-Pos-Pathogen_Box_20200109-processed.mgf',
        'tony_MSMS-Pos-PlaSMA-processed.mgf',
        'tony_MSMS-Pos-Respect-processed.mgf',
        'tony_MSMS-Pos-Vaniya-Fiehn_Natural_Products_Library_20200109-processed.mgf',
        'tony_MSMS-Pos-bmdms-np_20200811-processed.mgf',
        
        'tony_RIKEN_LIPIDOMICS_processed.mgf','tony_RIKEN_PlaSMA_processed.mgf','tony_RIKEN_RESPECT_processed.mgf',
        
        'MassBank_NIST.mgf', 'matchms_select.mgf', 'tony_FoxRamos2019-processed.mgf', 'MoNA-export-LC-MS-MS_Spectra.mgf', 'tony_MoNA_ShuoHan2021.mgf',
       ]

dup_df['INSTRUMENT_ALL'] = (
    dup_df['INSTRUMENT_TYPE'].fillna('') + ';' +
    dup_df['INSTRUMENT'].fillna('') + ';' +
    dup_df['SOURCE_INSTRUMENT'].fillna('') + ';' +
    (dup_df['INSTRUMENT_MANUFACTURER_AND_MODEL'].fillna('') if 'INSTRUMENT_MANUFACTURER_AND_MODEL' in dup_df.columns else '')
)

columns_to_replace = ['INSTRUMENT_ALL']
dup_df[columns_to_replace] = dup_df[columns_to_replace].replace(r'^\s*$', np.nan, regex=True)
dup_df[columns_to_replace] = dup_df[columns_to_replace].replace(r'^;+$', np.nan, regex=True)

dup_df['INSTRUMENT_MAIN'] = np.nan

instrument_mapping = {
    'qtof': 'qtof',
    'q-tof': 'qtof',
    'qtfo': 'qtof',
    'orbitrap': 'orbitrap',
    'qqq': 'qqq',
    'qq': 'qqq',
    'ittof': 'iontrap',
    'it-tof': 'iontrap',
    'ion trap': 'iontrap',
    'ion-trap': 'iontrap',
    'tripletof': 'qtof',
    'triple tof': 'qtof',
    'qft': 'orbitrap',
    'itft': 'orbitrap',
    'hf': 'orbitrap',
    'q exactive': 'orbitrap',
    'q-exactive': 'orbitrap',
    'ttof': 'qtof',
    'LC-ESI-CID; Velos': 'iontrap',
    'LC-ESI-HCD Velos': 'orbitrap',
    'Bruker timsTOF Pro': 'qtof',
    'LC-ESI-QIT;4000Q TRAP': 'qtof',
    'LC-ESI-TOF impact HD': 'qtof',
    'LC-ESI-HCD Lumos': 'orbitrap',
    'LC-ESI-CID; Lumos': 'iontrap',
    'FTMS-ESI': 'orbitrap',
    'Thermo Finnigan LTQ': 'iontrap',
    # 'Hybrid FT': 'iontrap/orbitrap',
    'MALDI-TOFTOF JMS-S3000': 'qtof',
    'ESI-TOF;micrOTOF-Q': 'qtof',
    'ESI-HCD': 'orbitrap',
    'impact HD': 'qtof',
    'MALDI-TOFTOF': 'qtof',
    'LC-ESI-TOF;LCT Micromass': 'qtof',
    'FAB-EBEB': 'qtof',
    'MALDI-QIT;AXIMA QIT': 'qtof',
    'QIT;API QSTAR': 'qtof',
    'Waters SYNAPT': 'qtof',
    'HCD': 'orbitrap'
}

# Apply the mapping step by step
for pattern, main_type in instrument_mapping.items():
    dup_df.loc[
        (dup_df['INSTRUMENT_ALL'].str.contains(pattern, case=False, na=False)) & 
        (dup_df['INSTRUMENT_MAIN'].isnull()), 
        'INSTRUMENT_MAIN'
    ] = main_type

dup_df.loc[(dup_df['COLLISION_ENERGY'].str.contains('HCD', na=False)) & (dup_df['INSTRUMENT_MAIN'].isna()), 'INSTRUMENT_MAIN'] = 'orbitrap'

instrument_order = [np.nan, 'orbitrap', 'qtof', 'iontrap', 'qqq']

def get_instrument_order(x):
    if pd.isna(x):
        return instrument_order.index(np.nan)
    # Otherwise use normal index lookup
    return instrument_order.index(x) if x in instrument_order else len(instrument_order)

dup_df['instrument_order'] = dup_df['INSTRUMENT_MAIN'].apply(get_instrument_order)

dup_df['source_order'] = dup_df['SOURCE'].apply(
    lambda x: library_order.index(x) if x in library_order else len(library_order)
)

dup_df = dup_df.sort_values(
    by=['instrument_order', 'source_order'], 
    ascending=[True, False]
).drop(columns=['instrument_order', 'source_order'])

dup_df['NUM_PEAKS'] = pd.to_numeric(dup_df['NUM_PEAKS'])

dup_df.drop(
    columns=[col for col in dup_df.columns if 'INSTRUMENT' in col or 'COLLISION' in col],
    inplace=True
)

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
