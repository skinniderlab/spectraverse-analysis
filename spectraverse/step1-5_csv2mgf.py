import pandas as pd
import numpy as np
import re, os, sys, math
from rdkit import Chem
from concurrent.futures import ProcessPoolExecutor, as_completed
import requests
from rdkit.Chem import Descriptors
from concurrent.futures import ThreadPoolExecutor, as_completed

print("Step1-5: CSV to MGF conversion (accomplished with metadata modification)")

in_csv_dir = sys.argv[1]
in_mgf_dir = sys.argv[2]

benchmark_dir = os.path.dirname(in_csv_dir)

out_csv_dir = sys.argv[3]
out_mgf_dir = sys.argv[4]

metadata = pd.read_csv(in_csv_dir)

def combine_columns(row):
    if pd.isna(row['INCHIKEY']) and not pd.isna(row['INCHIAUX']):
        return row['INCHIAUX']
    elif not pd.isna(row['INCHIKEY']) and pd.isna(row['INCHIAUX']):
        return row['INCHIKEY']
    elif not pd.isna(row['INCHIKEY']) and not pd.isna(row['INCHIAUX']):
        assert row['INCHIKEY'] == row['INCHIAUX'], f"Values in columns INCHIKEY and INCHIAUX are not equal: {row['INCHIKEY']} != {row['INCHIAUX']}"
        return row['INCHIKEY']
    else:
        return np.nan
    
metadata['INCHIKEY_FINAL'] = metadata.apply(combine_columns, axis=1)  

inchi_pattern = re.compile(r'^InChI=1S/')
inchikey_pattern = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')

def classify_and_move_inchi(row):
    content = str(row['INCHI'])
    if inchi_pattern.match(content):
        row['INCHI'] = content
    elif content.startswith('1S/'):
        row['INCHI'] = 'InChI=' + content
    elif inchikey_pattern.match(content):
        row['INCHIKEY_FINAL'] = content
        row['INCHI'] = None
    else:
        row['SMILES'] = content
        row['INCHI'] = None    
    return row

rows_to_check = metadata[(metadata['SMILES'].isna())&(metadata['INCHI'].notna())]
metadata.loc[rows_to_check.index] = rows_to_check.apply(classify_and_move_inchi, axis=1)

def remove_dot_suffix(s):
    if isinstance(s, str):
        while s.endswith('.'):
            s = s[:-1]
    return s

def remove_nextline_suffix(s):
    if isinstance(s, str):
        while s.endswith('&gt;'):
            s = s[:-4]
    return s   

def remove_semicolon_suffix(s):
    if isinstance(s, str):
        s = s.rstrip(';') 
        s = s.lstrip(';') 
    return s

def remove_between_pipes(s):
    if isinstance(s, str):
        return re.sub(r'\|.*?\|', '', s)
    return s

def add_inchi_prefix(s):
    if not isinstance(s, str):
        return s

    inchi = s
    if s.startswith('1S/'):
        inchi = 'InChI=' + s
    elif inchikey_pattern.match(s):
        inchi = s

    if inchi_pattern.match(inchi):
        try:
            mol = Chem.MolFromInchi(inchi)
            if mol:
                return Chem.MolToSmiles(mol)
        except Exception:
            pass
        return inchi
    return s

metadata['SMILES'] = metadata['SMILES'].str.strip().str.replace(' ', '')
metadata['SMILES'] = metadata['SMILES'].apply(remove_dot_suffix)
metadata['SMILES'] = metadata['SMILES'].apply(remove_nextline_suffix)
metadata['SMILES'] = metadata['SMILES'].apply(remove_semicolon_suffix)
metadata['SMILES'] = metadata['SMILES'].apply(remove_between_pipes)
metadata['SMILES'] = metadata['SMILES'].apply(add_inchi_prefix)

metadata['INCHIKEY_FINAL'] = metadata['INCHIKEY_FINAL'].replace('N/A', np.nan)

smiles_inchikey_path = benchmark_dir + '/inchikey_smiles.csv'
smiles_inchikey_df = pd.read_csv(smiles_inchikey_path)

metadata = metadata.merge(smiles_inchikey_df, left_on='INCHIKEY_FINAL', right_on='INCHIKEY', how='left')
metadata.rename(columns={'SMILES_x': 'SMILES', 'INCHIKEY_x':'INCHIKEY'}, inplace=True)
rows_to_check = metadata[(metadata['SMILES'].isna())&(metadata['INCHI'].isna())&(metadata['INCHIKEY_FINAL'].notna())]
metadata.loc[rows_to_check.index, 'SMILES'] = metadata.loc[rows_to_check.index, 'SMILES_y']
metadata = metadata.drop(columns=['SMILES_y', 'INCHIKEY_y'])

metadata[metadata['SMILES']=='H][C@](O)(CO)COP(O)(=O)OC[C@@]([H])(COC(=O)CCCCCCC\C=C/CCCCCC)OC(=O)CCCCCCC\C=C/CCCCCCCC'] == '[H][C@](O)(CO)COP(O)(=O)OC[C@@]([H])(COC(=O)CCCCCCC\C=C/CCCCCC)OC(=O)CCCCCCC\C=C/CCCCCCCC'

metadata['PRECURSOR_MZ'] = metadata['PRECURSOR_MZ'].fillna(-1)

smiles_na_cleaned_file = benchmark_dir + '/no_smiles_link_new.csv'
ref = pd.read_csv(smiles_na_cleaned_file)

smiles_na_unique = ref['COMPOUND_NAME'].unique()
for i, compound in enumerate(smiles_na_unique):
    metadata.loc[(metadata['SMILES'].isna()) & (metadata['COMPOUND_NAME'] == compound), 'SMILES'] = ref['SMILES'][i]

def get_compound_name_from_gnps(spectrum_id):
    print(f"Processing spectrum ID {spectrum_id}")
    url = f"https://gnps.ucsd.edu/ProteoSAFe/SpectrumCommentServlet?SpectrumID={spectrum_id}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'annotations' in data and len(data['annotations']) > 0:
            compound_name = data['annotations'][0].get('Compound_Name', None)
            return compound_name
        else:
            print(f"No annotations found for spectrum ID {spectrum_id}")
            return None
    else:
        print(f"Failed to retrieve data for spectrum ID {spectrum_id}")
        return None

def process_row(row):
    if row['SOURCE'] == 'matchms_select.mgf' and pd.isna(row['SMILES']) and (pd.isna(row['INCHI'])):
        return get_compound_name_from_gnps(row['SPECTRUMID'])
    else:
        return row['COMPOUND_NAME']

with ThreadPoolExecutor() as executor:
    futures = {executor.submit(process_row, row): index for index, row in metadata.iterrows()}
    for future in as_completed(futures):
        index = futures[future]
        try:
            result = future.result()
            metadata.at[index, 'COMPOUND_NAME'] = result
        except Exception as e:
            print(f"Error processing row {index}: {e}")

def calculate_exact_mass(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        exact_mass = Descriptors.ExactMolWt(mol)
    except:
        exact_mass = None    
    return exact_mass

metadata['PARENT_MASS'] = metadata['SMILES'].apply(calculate_exact_mass)

metadata['MASS_DIFF'] =  metadata['PRECURSOR_MZ'] - metadata['PARENT_MASS']

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

        smiles = metadata['SMILES'][i]
        if pd.isna(smiles):
            mgf_file.write("SMILES=\n")
        else:
            mgf_file.write("SMILES={}\n".format(smiles))

        inchi = metadata['INCHI'][i]
        if pd.isna(inchi):
            mgf_file.write("INCHI=\n")
        else:
            mgf_file.write("INCHI={}\n".format(inchi))    

        inchikey = metadata['INCHIKEY_FINAL'][i]
        if pd.isna(inchikey):
            mgf_file.write("INCHIKEY=\n")
        else:
            mgf_file.write("INCHIKEY={}\n".format(inchikey))  

        mgf_file.write("IONMODE={}\n".format(metadata['IONMODE'][i]))     

        adduct = metadata['ADDUCT'][i]
        if pd.isna(adduct):
            mgf_file.write("ADDUCT=\n")
        else:
            mgf_file.write("ADDUCT={}\n".format(adduct))

        # compound_name = metadata['COMPOUND_NAME'][i]
        # if not pd.isna(compound_name):
        #     mgf_file.write("COMPOUND_NAME={}\n".format(compound_name))    
            
        mgf_file.write("NUM_PEAKS={}\n".format(int(metadata['NUM_PEAKS'][i])))
        mgf_file.write("PRECURSOR_MZ={}\n".format(metadata['PRECURSOR_MZ'][i]))

        mgf_file.write("MS_LEVEL={}\n".format(metadata['MS_LEVEL'][i]))

        for line in info[i]:
            if re.match(r'^\d', line):
                mgf_file.write(line+'\n')

        mgf_file.write("END IONS\n\n")        

metadata.to_csv(out_csv_dir, index=False, compression='gzip')
