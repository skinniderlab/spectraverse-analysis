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

print("Step3-7: Adding INSTRUMENT and COLLISION_ENERGY information back to the processed mgf file and finalizing the metadata")

mgf_orig = sys.argv[1]
mgf_processed = sys.argv[2]
 
output_csv_dir = sys.argv[3]
output_mgf_dir = sys.argv[4]


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

info_orig = extract_info(mgf_orig)
info_processed = extract_info(mgf_processed)

info_orig_dict = {
    (next(item for item in sublist if item.startswith('INDEX=')),
     next(item for item in sublist if item.startswith('SOURCE='))): sublist
    for sublist in info_orig
}

for sublist2 in info_processed:
    index2 = [item for item in sublist2 if item.startswith('INDEX=')][0]
    source2 = [item for item in sublist2 if item.startswith('SOURCE=')][0]

    matching_sublist1 = info_orig_dict.get((index2, source2))
    if matching_sublist1 is None:
        print(f"No match found for index2: {index2}, source2: {source2}")
        print(f"sublist2: {sublist2}")

    assert matching_sublist1 is not None, f"No match found for index2: {index2}, source2: {source2}"
    
    for item in matching_sublist1:
        if 'INSTRUMENT' in item or 'COLLISION' in item:
            sublist2.insert(1, item)  

benchmark_dir = os.path.dirname(mgf_processed)

temp_mgf = benchmark_dir  + '/test-can-rem-lowres-prec-highfragmass-matchms-dotcheck-noiserem-temp.mgf'

with open(temp_mgf, 'w') as file:
    for sublist in info_processed:
        for item in sublist:
            file.write("%s\n" % item)
        file.write("\n")


rows = []  # Initialize an empty list to hold the rows

# Initialize an empty dictionary to hold the current row
row = {}

# Open the file and read it line by line
with open(temp_mgf, 'r') as file:
    for line in file:
        # Strip the newline character from the end of the line
        line = line.rstrip('\n')

        # If the line is "END IONS", add the current row to the list and reset the row
        if line == 'END IONS':
            rows.append(row)
            row = {}
        # If the line contains an equals sign and the first part of the line is not a number, split it into a column name and a value and add them to the row
        elif '=' in line and not line.split('=')[0].strip().isdigit():
            column, value = line.split('=', 1)
            row[column] = value

df = pd.DataFrame(rows)

def process_charge(row):
    if row['CHARGE'] == '0':
        if row['IONMODE'] == 'positive':
            return 1
        elif row['IONMODE'] == 'negative':
            return -1
    else:
        if row['CHARGE'] == '1-':
            return -1
        elif row['CHARGE'] == '1+':
            return 1

# Apply the function to the DataFrame
df['CHARGE'] = df.apply(process_charge, axis=1)

instrument_columns = [
    'INSTRUMENT_TYPE', 'INSTRUMENT', 'SOURCE_INSTRUMENT', 'INSTRUMENT_MANUFACTURER_AND_MODEL'
]
available_instrument_columns = [col for col in instrument_columns if col in df.columns]

if available_instrument_columns:
    df['INSTRUMENT_ALL'] = df[available_instrument_columns].fillna('').agg(';'.join, axis=1)
else:
    print("No INSTRUMENT columns are available. Skipping INSTRUMENT_ALL.")

collision_energy_columns = [
    'COLLISION_ENERGY', 'COLLISION_ENERGY_3', 'COLLISION_ENERGY_2',
    'COLLISION_ENERGY_VOLTAGE', 'NOMRALIZED_COLLISION_ENERGY',
    'NORMALIZED_COLLISION_ENERGY', 'COLLISION_ENERGY_SPREAD'
]
available_collision_energy_columns = [col for col in collision_energy_columns if col in df.columns]

if available_collision_energy_columns:
    df['COLLISION_ENERGY_ALL'] = df[available_collision_energy_columns].fillna('').agg(';'.join, axis=1)
else:
    print("No COLLISION_ENERGY columns are available. Skipping COLLISION_ENERGY_ALL.")


columns_to_replace = ['INSTRUMENT_ALL', 'COLLISION_ENERGY_ALL']
df[columns_to_replace] = df[columns_to_replace].replace(r'^\s*$', np.nan, regex=True)
df[columns_to_replace] = df[columns_to_replace].replace(r'^;+$', np.nan, regex=True)


df['INSTRUMENT_MAIN'] = np.nan

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

for pattern, main_type in instrument_mapping.items():
    df.loc[
        (df['INSTRUMENT_ALL'].str.contains(pattern, case=False, na=False)) & 
        (df['INSTRUMENT_MAIN'].isnull()), 
        'INSTRUMENT_MAIN'
    ] = main_type

df.loc[(df['COLLISION_ENERGY'].str.contains('HCD', na=False)) & (df['INSTRUMENT_MAIN'].isna()), 'INSTRUMENT_MAIN'] = 'orbitrap'

qqq_df = df[df['INSTRUMENT_MAIN'] == 'qqq']
df = df[df['INSTRUMENT_MAIN'] != 'qqq']

df['COLLISION_ENERGY_MAIN1'] = np.nan
df['COLLISION_ENERGY_MAIN2'] = np.nan
df['COLLISION_ENERGY_MAIN3'] = np.nan
df['N_COLLISION_ENERGY_MAIN1'] = np.nan
df['N_COLLISION_ENERGY_MAIN2'] = np.nan
df['N_COLLISION_ENERGY_MAIN3'] = np.nan

def process_collision_energy_ce1(row):
    if isinstance(row['COLLISION_ENERGY'], str) and re.match(r'^\d+(\.\d+)?$', row['COLLISION_ENERGY']):
        if row['INSTRUMENT_MAIN'] in ['orbitrap', 'iontrap']:
            row['N_COLLISION_ENERGY_MAIN1'] = float(row['COLLISION_ENERGY'])
        else:
            row['COLLISION_ENERGY_MAIN1'] = float(row['COLLISION_ENERGY'])
    elif isinstance(row['COLLISION_ENERGY'], str) and re.match(r'^\s*(\d+(\.\d+)?)\s*(eV|V|ev)?\s*$', row['COLLISION_ENERGY'], re.IGNORECASE) and all(keyword not in row['COLLISION_ENERGY'] for keyword in ['NCE', 'HCD', 'CID']):
        row['COLLISION_ENERGY_MAIN1'] = float(re.search(r'\d+(\.\d+)?', row['COLLISION_ENERGY']).group())
    return row

df = df.apply(process_collision_energy_ce1, axis=1)

df.loc[
    (df['COLLISION_ENERGY'] == '1030') & (df['INSTRUMENT_MAIN'].isin(['orbitrap', 'iontrap'])),
    ['N_COLLISION_ENERGY_MAIN1', 'N_COLLISION_ENERGY_MAIN2']
] = [10, 30]

df.loc[
    (df['COLLISION_ENERGY'] == '1030') & (~df['INSTRUMENT_MAIN'].isin(['orbitrap', 'iontrap'])),
    ['COLLISION_ENERGY_MAIN1', 'COLLISION_ENERGY_MAIN2']
] = [10, 30]

df.loc[
    (df['COLLISION_ENERGY'] == '10, 30') & (df['INSTRUMENT_MAIN'].isin(['orbitrap', 'iontrap'])),
    ['N_COLLISION_ENERGY_MAIN1', 'N_COLLISION_ENERGY_MAIN2']
] = [10, 30]

df.loc[
    (df['COLLISION_ENERGY'] == '10, 30') & (~df['INSTRUMENT_MAIN'].isin(['orbitrap', 'iontrap'])),
    ['COLLISION_ENERGY_MAIN1', 'COLLISION_ENERGY_MAIN2']
] = [10, 30]

df.loc[
    (df['COLLISION_ENERGY'] == '25,60,100 (stepped)') & (df['INSTRUMENT_MAIN'].isin(['orbitrap', 'iontrap'])),
    ['N_COLLISION_ENERGY_MAIN1', 'N_COLLISION_ENERGY_MAIN2', 'N_COLLISION_ENERGY_MAIN3']
] = [25, 60, 100]

df.loc[
    (df['COLLISION_ENERGY'] == '25,60,100 (stepped)') & (~df['INSTRUMENT_MAIN'].isin(['orbitrap', 'iontrap'])),
    ['COLLISION_ENERGY_MAIN1', 'COLLISION_ENERGY_MAIN2', 'COLLISION_ENERGY_MAIN3']
] = [25, 60, 100]

df.loc[
    (df['COLLISION_ENERGY'] == 'hcd30.00') & (df['INSTRUMENT_MAIN'].isin(['orbitrap', 'iontrap'])),
    ['N_COLLISION_ENERGY_MAIN1']
] = [30]

df.loc[
    (df['COLLISION_ENERGY'] == 'hcd30.00') & (~df['INSTRUMENT_MAIN'].isin(['orbitrap', 'iontrap'])),
    ['COLLISION_ENERGY_MAIN1']
] = [30]

df.loc[df['COLLISION_ENERGY'] == '-35 V', ['COLLISION_ENERGY_MAIN1']] = [35]
df.loc[df['COLLISION_ENERGY'] == '35 eV FT-MS II', ['COLLISION_ENERGY_MAIN1']] = [35]
df.loc[df['COLLISION_ENERGY'] == '35 eV IT-MS', ['COLLISION_ENERGY_MAIN1']] = [35]
df.loc[df['COLLISION_ENERGY'] == '65 eV FT-MS', ['COLLISION_ENERGY_MAIN1']] = [65]
df.loc[df['COLLISION_ENERGY'] == '35 eV FT-MS', ['COLLISION_ENERGY_MAIN1']] = [35]
df.loc[df['COLLISION_ENERGY'] == '45 eV (ion-source)', ['COLLISION_ENERGY_MAIN1']] = [45]
df.loc[df['COLLISION_ENERGY'] == '30 V (ramped)', ['COLLISION_ENERGY_MAIN1']] = [30]

df.loc[df['NORMALIZED_COLLISION_ENERGY'] == '40', ['N_COLLISION_ENERGY_MAIN1']] = [40]
df.loc[df['NOMRALIZED_COLLISION_ENERGY'] == '40', ['N_COLLISION_ENERGY_MAIN1']] = [40]
df.loc[df['NORMALIZED_COLLISION_ENERGY'] == '10,20,40', ['N_COLLISION_ENERGY_MAIN1', 'N_COLLISION_ENERGY_MAIN2', 'N_COLLISION_ENERGY_MAIN3']] = [10, 20, 40]
df.loc[df['NOMRALIZED_COLLISION_ENERGY'] == '10,20,40', ['N_COLLISION_ENERGY_MAIN1', 'N_COLLISION_ENERGY_MAIN2', 'N_COLLISION_ENERGY_MAIN3']] = [10, 20, 40]



def process_collision_energy_ce2(row):
    if isinstance(row['COLLISION_ENERGY_2'], str) and re.match(r'^\d+(\.\d+)?$', row['COLLISION_ENERGY_2']):
        if row['INSTRUMENT_MAIN'] in ['orbitrap', 'iontrap']:
            row['N_COLLISION_ENERGY_MAIN2'] = float(row['COLLISION_ENERGY_2'])
        else:
            row['COLLISION_ENERGY_MAIN2'] = float(row['COLLISION_ENERGY_2'])
    elif isinstance(row['COLLISION_ENERGY_2'], str) and re.match(r'^\s*(\d+(\.\d+)?)\s*(eV|V|ev|HCD|CID)?\s*$', row['COLLISION_ENERGY_2'], re.IGNORECASE) and 'NCE' not in row['COLLISION_ENERGY_2']:
        row['COLLISION_ENERGY_MAIN2'] = float(re.search(r'\d+(\.\d+)?', row['COLLISION_ENERGY_2']).group())
    return row

if 'COLLISION_ENERGY_2' not in df.columns:
    print("Column 'COLLISION_ENERGY_2' not found in the DataFrame.")
else:
    df = df.apply(process_collision_energy_ce2, axis=1)

def process_collision_energy_ce3(row):
    if isinstance(row['COLLISION_ENERGY_3'], str) and re.match(r'^\d+(\.\d+)?$', row['COLLISION_ENERGY_3']):
        if row['INSTRUMENT_MAIN'] in ['orbitrap', 'iontrap']:
            row['N_COLLISION_ENERGY_MAIN3'] = float(row['COLLISION_ENERGY_3'])
        else:
            row['COLLISION_ENERGY_MAIN3'] = float(row['COLLISION_ENERGY_3'])
    elif isinstance(row['COLLISION_ENERGY_3'], str) and re.match(r'^\s*(\d+(\.\d+)?)\s*(eV|V|ev|HCD|CID)?\s*$', row['COLLISION_ENERGY_3'], re.IGNORECASE) and 'NCE' not in row['COLLISION_ENERGY_3']:
        row['COLLISION_ENERGY_MAIN3'] = float(re.search(r'\d+(\.\d+)?', row['COLLISION_ENERGY_3']).group())
    return row

if 'COLLISION_ENERGY_3' not in df.columns:
    print("Column 'COLLISION_ENERGY_3' not found in the DataFrame.")
else:
    df = df.apply(process_collision_energy_ce3, axis=1)

def process_other_columns_ce(df, column_name):
    mask = df[column_name].notna()
    processed = df.loc[mask, column_name].str.replace(r'[eE]?[vV]', '', regex=True).str.split(',')
    expanded = processed.apply(lambda x: pd.Series([float(v.strip()) for v in x] + [np.nan] * (3 - len(x))))
    df.loc[mask, ['COLLISION_ENERGY_MAIN1', 'COLLISION_ENERGY_MAIN2', 'COLLISION_ENERGY_MAIN3']] = expanded.values
    return df

if 'COLLISION_ENERGY_VOLTAGE' not in df.columns:
    print("Column 'COLLISION_ENERGY_VOLTAGE' not found in the DataFrame.")
else:
    df = process_other_columns_ce(df, 'COLLISION_ENERGY_VOLTAGE')

    df.loc[df['COLLISION_ENERGY_VOLTAGE'] == '30', ['COLLISION_ENERGY_MAIN1']] = np.nan
    df.loc[df['COLLISION_ENERGY_VOLTAGE'] == '30', ['N_COLLISION_ENERGY_MAIN1']] = [30]

def extract_ramp_values(value, original_main1, original_main2):
    if isinstance(value, str) and re.search(r'ramp', value, re.IGNORECASE):
        match = re.search(r'(\d+(\.\d+)?)\s*-\s*(\d+(\.\d+)?)\s*(eV|V)', value, re.IGNORECASE)
        if match:
            return round(float(match.group(1)), 6), round(float(match.group(3)), 6)
    return original_main1, original_main2

df[['COLLISION_ENERGY_MAIN1', 'COLLISION_ENERGY_MAIN2']] = df.apply(
    lambda row: pd.Series(
        extract_ramp_values(row['COLLISION_ENERGY'], row['COLLISION_ENERGY_MAIN1'], row['COLLISION_ENERGY_MAIN2'])
    ),
    axis=1
)

def process_hcd_txt(row):
    if isinstance(row['COLLISION_ENERGY'], str) and re.match(r'^\d+\s*(?:HCD|CID)$', row['COLLISION_ENERGY'], re.IGNORECASE):
        output = re.search(r'\d+', row['COLLISION_ENERGY']).group()
        result = round(float(output), 6)
        if row['INSTRUMENT_MAIN'] in ['orbitrap', 'iontrap']:
            row['N_COLLISION_ENERGY_MAIN1'] = result
        else:
            row['COLLISION_ENERGY_MAIN1'] = result
    return row

df = df.apply(process_hcd_txt, axis=1)

def process_ce_txt(row):
    if isinstance(row['COLLISION_ENERGY'], str) and re.match(r'^CE\s*\d+$', row['COLLISION_ENERGY'], re.IGNORECASE):
        output = re.search(r'\d+', row['COLLISION_ENERGY']).group()
        return round(float(output), 6)
    else:
        return row['COLLISION_ENERGY_MAIN1']
        
df['COLLISION_ENERGY_MAIN1'] = df.apply(process_ce_txt, axis=1)

def extract_range_values(value, original_main1, original_main2):
    if isinstance(value, str):
        match = re.match(r'^\s*(\d+(\.\d+)?)\s*[-|>]+\s*(\d+(\.\d+)?)\s*(eV|V)?\s*$', value, re.IGNORECASE)
        if match:
            return round(float(match.group(1)), 6), round(float(match.group(3)), 6)
    return original_main1, original_main2

df[['COLLISION_ENERGY_MAIN1', 'COLLISION_ENERGY_MAIN2']] = df.apply(
    lambda row: pd.Series(
        extract_range_values(row['COLLISION_ENERGY'], row['COLLISION_ENERGY_MAIN1'], row['COLLISION_ENERGY_MAIN2'])
    ),
    axis=1
)

def extract_float_from_string(value):
    match = re.search(r'\d+(\.\d+)?', value)
    return float(match.group())

def calculate_ce(nce, precursor_mz):
    if isinstance(precursor_mz, str):
        precursor_mz = extract_float_from_string(precursor_mz)
    else:
        precursor_mz = float(precursor_mz)

    return (nce * precursor_mz)  / 500

def process_single_nce(row):
    value = row['COLLISION_ENERGY']
    precursor_mz = row['PRECURSOR_MZ']
    
    if pd.isnull(row['COLLISION_ENERGY_MAIN1']) and isinstance(value, str):
        # Exclude invalid cases
        if any(keyword in value.lower() for keyword in ['stepped', 'ramp', 'by', ',', 'scaled', 'kv', 'argon', '-', 'pqd']):
            return row['COLLISION_ENERGY_MAIN1'], row['N_COLLISION_ENERGY_MAIN1']
        
        match = re.match(r'.*?(\d+(\.\d+)?)\s*(\(|)?\s*(%|NCE|nominal)?\s*(\)|)?.*', value, re.IGNORECASE)
        if match:
            nce = float(match.group(1))
            return round(calculate_ce(nce, precursor_mz), 6), round(nce, 6)
        
        if value == "70(":
            nce = 70
            return round(calculate_ce(nce, precursor_mz), 6), round(nce, 6)
    
    return row['COLLISION_ENERGY_MAIN1'], row['N_COLLISION_ENERGY_MAIN1']

df[['COLLISION_ENERGY_MAIN1', 'N_COLLISION_ENERGY_MAIN1']] = df.apply(process_single_nce, axis=1, result_type='expand')

def process_stepped_nce(row):
    value = row['COLLISION_ENERGY']
    precursor_mz = row['PRECURSOR_MZ'] 
    
    if pd.isnull(row['COLLISION_ENERGY_MAIN1']) and isinstance(value, str) and 'NCE' in value:
        match = re.search(r'(\d+(\.\d+)?)\s*[-,]\s*(\d+(\.\d+)?)\s*[-,]\s*(\d+(\.\d+)?)', value)
        if match:
            nce_num1, nce_num2, nce_num3 = float(match.group(1)), float(match.group(3)), float(match.group(5))
            ce_num1 = calculate_ce(nce_num1, precursor_mz)
            ce_num2 = calculate_ce(nce_num2, precursor_mz)
            ce_num3 = calculate_ce(nce_num3, precursor_mz)
            return pd.Series([round(ce_num1, 6), round(ce_num2, 6), round(ce_num3, 6), round(nce_num1, 6), round(nce_num2, 6), round(nce_num3, 6)])
    
    return pd.Series([
         row['COLLISION_ENERGY_MAIN1'], 
         row['COLLISION_ENERGY_MAIN2'], 
         row['COLLISION_ENERGY_MAIN3'],
         row['N_COLLISION_ENERGY_MAIN1'],
         row['N_COLLISION_ENERGY_MAIN2'],
         row['N_COLLISION_ENERGY_MAIN3']
         ])

df[['COLLISION_ENERGY_MAIN1', 
    'COLLISION_ENERGY_MAIN2', 
    'COLLISION_ENERGY_MAIN3',
    'N_COLLISION_ENERGY_MAIN1',
    'N_COLLISION_ENERGY_MAIN2',
    'N_COLLISION_ENERGY_MAIN3']] = df.apply(process_stepped_nce, axis=1)

for index, precursor_mz in df.loc[df['COLLISION_ENERGY'] == 'Ramp 20%-70% (nominal)', 'PRECURSOR_MZ'].items():
    ce1 = calculate_ce(20, precursor_mz)
    ce2 = calculate_ce(70, precursor_mz)
    df.loc[index, ['COLLISION_ENERGY_MAIN1', 'COLLISION_ENERGY_MAIN2', 'N_COLLISION_ENERGY_MAIN1', 'N_COLLISION_ENERGY_MAIN2']] = [ce1, ce2, 20, 70]

for index, precursor_mz in df.loc[df['COLLISION_ENERGY'] == '42% (nominal) with stepped collision energies 30% and 55%', 'PRECURSOR_MZ'].items():
    ce1 = calculate_ce(30, precursor_mz)
    ce2 = calculate_ce(42, precursor_mz)
    ce3 = calculate_ce(55, precursor_mz)
    df.loc[index, ['COLLISION_ENERGY_MAIN1', 'COLLISION_ENERGY_MAIN2', 'COLLISION_ENERGY_MAIN3', 'N_COLLISION_ENERGY_MAIN1', 'N_COLLISION_ENERGY_MAIN2', 'N_COLLISION_ENERGY_MAIN3']] = [ce1, ce2, ce3, 30, 42, 55]

for index, precursor_mz in df.loc[df['COLLISION_ENERGY'] == '25,60,100% (stepped)', 'PRECURSOR_MZ'].items():
    ce1 = calculate_ce(25, precursor_mz)
    ce2 = calculate_ce(60, precursor_mz)
    ce3 = calculate_ce(100, precursor_mz)
    df.loc[index, ['COLLISION_ENERGY_MAIN1', 'COLLISION_ENERGY_MAIN2', 'COLLISION_ENERGY_MAIN3', 'N_COLLISION_ENERGY_MAIN1', 'N_COLLISION_ENERGY_MAIN2', 'N_COLLISION_ENERGY_MAIN3']] = [ce1, ce2, ce3, 25, 60, 100]

for index, precursor_mz in df.loc[df['COLLISION_ENERGY'] == '20,40,130% (stepped)', 'PRECURSOR_MZ'].items():
    ce1 = calculate_ce(20, precursor_mz)
    ce2 = calculate_ce(40, precursor_mz)
    ce3 = calculate_ce(130, precursor_mz)
    df.loc[index, ['COLLISION_ENERGY_MAIN1', 'COLLISION_ENERGY_MAIN2', 'COLLISION_ENERGY_MAIN3', 'N_COLLISION_ENERGY_MAIN1', 'N_COLLISION_ENERGY_MAIN2', 'N_COLLISION_ENERGY_MAIN3']] = [ce1, ce2, ce3, 20, 40, 130]

for index, precursor_mz in df.loc[df['COLLISION_ENERGY'] == 'HCD stepped 10, 30, 50 eV', 'PRECURSOR_MZ'].items():
    ce1 = calculate_ce(10, precursor_mz)
    ce2 = calculate_ce(30, precursor_mz)
    ce3 = calculate_ce(50, precursor_mz)
    df.loc[index, ['COLLISION_ENERGY_MAIN1', 'COLLISION_ENERGY_MAIN2', 'COLLISION_ENERGY_MAIN3', 'N_COLLISION_ENERGY_MAIN1', 'N_COLLISION_ENERGY_MAIN2', 'N_COLLISION_ENERGY_MAIN3']] = [ce1, ce2, ce3, 10, 30, 50]


def calculate_nce(nce, precursor_mz):
    if isinstance(precursor_mz, str):
        precursor_mz = extract_float_from_string(precursor_mz)
    else:
        precursor_mz = float(precursor_mz)

    return (nce * 500) / precursor_mz

def process_missing_nce(row, input_column, output_column, precursor_column):
    value = row[input_column]
    precursor_mz = row[precursor_column]
    
    if pd.notnull(value) and pd.isnull(row[output_column]) and isinstance(value, float):
        ce = float(value)
        nce = calculate_nce(ce, precursor_mz)
        return round(nce, 6)

    return row[output_column]

df['N_COLLISION_ENERGY_MAIN1'] = df.apply(process_missing_nce, axis=1, args=('COLLISION_ENERGY_MAIN1', 'N_COLLISION_ENERGY_MAIN1', 'PRECURSOR_MZ'))  
df['N_COLLISION_ENERGY_MAIN2'] = df.apply(process_missing_nce, axis=1, args=('COLLISION_ENERGY_MAIN2', 'N_COLLISION_ENERGY_MAIN2', 'PRECURSOR_MZ'))  
df['N_COLLISION_ENERGY_MAIN3'] = df.apply(process_missing_nce, axis=1, args=('COLLISION_ENERGY_MAIN3', 'N_COLLISION_ENERGY_MAIN3', 'PRECURSOR_MZ'))  

def process_missing_ce(row, input_column, output_column, precursor_column):
    value = row[input_column]
    precursor_mz = row[precursor_column]
    
    if pd.notnull(value) and pd.isnull(row[output_column]) and isinstance(value, float):
        ce = float(value)
        nce = calculate_ce(ce, precursor_mz)
        return round(nce, 6)

    return row[output_column]

df['COLLISION_ENERGY_MAIN1'] = df.apply(process_missing_ce, axis=1, args=('N_COLLISION_ENERGY_MAIN1', 'COLLISION_ENERGY_MAIN1', 'PRECURSOR_MZ'))  
df['COLLISION_ENERGY_MAIN2'] = df.apply(process_missing_ce, axis=1, args=('N_COLLISION_ENERGY_MAIN2', 'COLLISION_ENERGY_MAIN2', 'PRECURSOR_MZ'))  
df['COLLISION_ENERGY_MAIN3'] = df.apply(process_missing_ce, axis=1, args=('N_COLLISION_ENERGY_MAIN3', 'COLLISION_ENERGY_MAIN3', 'PRECURSOR_MZ'))  

instrument_columns_to_drop = ['INSTRUMENT_TYPE','INSTRUMENT', 'SOURCE_INSTRUMENT', 'INSTRUMENT_MANUFACTURER_AND_MODEL', 'INSTRUMENT_ALL']

df = df.drop(columns=instrument_columns_to_drop, errors='ignore')         

ce_columns_to_drop = ['COLLISION_ENERGY', 'COLLISION_ENERGY_3', 'COLLISION_ENERGY_2', 'COLLISION_ENERGY_VOLTAGE', 
                      'NOMRALIZED_COLLISION_ENERGY', 'NORMALIZED_COLLISION_ENERGY', 'COLLISION_ENERGY_SPREAD',
                      'COLLISION_ENERGY_ALL', 'COLLISION_CELL_RF', 'COLLISION_GAS', 'COLLISION_CELL_GAS']

df = df.drop(columns=ce_columns_to_drop, errors='ignore')          

additional_columns_to_drop = ['NUM_PEAKS_OLD']
df = df.drop(columns=additional_columns_to_drop)

df = df.rename(columns={
    'INSTRUMENT_MAIN': 'INSTRUMENT_TYPE',
    'COLLISION_ENERGY_MAIN1': 'COLLISION_ENERGY_1',
    'COLLISION_ENERGY_MAIN2': 'COLLISION_ENERGY_2',
    'COLLISION_ENERGY_MAIN3': 'COLLISION_ENERGY_3',
    'N_COLLISION_ENERGY_MAIN1': 'NORMALIZED_COLLISION_ENERGY_1',
    'N_COLLISION_ENERGY_MAIN2': 'NORMALIZED_COLLISION_ENERGY_2',
    'N_COLLISION_ENERGY_MAIN3': 'NORMALIZED_COLLISION_ENERGY_3'
})


def extract_info_index(filename, indices):
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

df_index = df.index

temp_info = extract_info_index(temp_mgf, df_index)

df = df.reset_index(drop=True)
df.insert(0, 'TITLE', df.index.to_series().apply(lambda x: f"SPECTRAVERSE{x + 1:08d}"))
df = df.drop(columns=['SOURCE', 'INDEX'])


df.to_csv(output_csv_dir, index=False, encoding='ascii', errors='ignore')

df = pd.read_csv(output_csv_dir)


with open(output_mgf_dir, 'w') as mgf_file:
    for i in range(len(df)):
        mgf_file.write("BEGIN IONS\n")

        mgf_file.write("TITLE={}\n".format(df['TITLE'][i]))
        mgf_file.write("FORMULA={}\n".format(df['FORMULA'][i]))
        mgf_file.write("SMILES={}\n".format(df['SMILES'][i]))
        mgf_file.write("INCHI={}\n".format(df['INCHI'][i]))
        mgf_file.write("INCHIKEY={}\n".format(df['INCHIKEY'][i]))
        mgf_file.write("IONMODE={}\n".format(df['IONMODE'][i]))
        mgf_file.write("ADDUCT={}\n".format(df['ADDUCT'][i]))
        # mgf_file.write("COMPOUND_NAME={}\n".format(df['COMPOUND_NAME'][i]))  
        mgf_file.write("NUM_PEAKS={}\n".format(df['NUM_PEAKS'][i]))
        mgf_file.write("PRECURSOR_MZ={}\n".format(df['PRECURSOR_MZ'][i]))
        mgf_file.write("PEPMASS={}\n".format(df['PRECURSOR_MZ'][i]))
        mgf_file.write("PARENT_MASS={}\n".format(df['PARENT_MASS'][i]))
        mgf_file.write("MS_LEVEL={}\n".format(df['MS_LEVEL'][i]))
        mgf_file.write("CHARGE={}\n".format(df['CHARGE'][i]))
        mgf_file.write("INSTRUMENT_TYPE={}\n".format(df['INSTRUMENT_TYPE'][i]))
        mgf_file.write("COLLISION_ENERGY_1={}\n".format(df['COLLISION_ENERGY_1'][i]))
        mgf_file.write("COLLISION_ENERGY_2={}\n".format(df['COLLISION_ENERGY_2'][i]))
        mgf_file.write("COLLISION_ENERGY_3={}\n".format(df['COLLISION_ENERGY_3'][i]))
        mgf_file.write("NORMALIZED_COLLISION_ENERGY_1={}\n".format(df['NORMALIZED_COLLISION_ENERGY_1'][i]))
        mgf_file.write("NORMALIZED_COLLISION_ENERGY_2={}\n".format(df['NORMALIZED_COLLISION_ENERGY_2'][i]))
        mgf_file.write("NORMALIZED_COLLISION_ENERGY_3={}\n".format(df['NORMALIZED_COLLISION_ENERGY_3'][i]))

        for line in temp_info[i]:
            if re.match(r'^\d', line):
                mgf_file.write(line+'\n')

        mgf_file.write("END IONS\n\n")   
