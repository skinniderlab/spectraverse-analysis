import numpy as np
import os, glob, re
import pandas as pd
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdMolDescriptors
from tqdm import tqdm
from multiprocessing import Pool
import warnings

print("Step3-2: Removal of unwanted adducts and correction of [M]- to [M-H]-")

warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

input_csv_dir = sys.argv[1]
input_mgf_dir = sys.argv[2]

output_csv_dir = sys.argv[3]
output_mgf_dir = sys.argv[4]

removed_index = []

metadata = pd.read_csv(input_csv_dir)

remove_neutralloss_adduct = ['[M-H2O+Na]+','[M-2H2O+Na]+', '[M-3H2O+Na]+', '[M-H2O+NH4]+', '[M+H-3H2O]+', '[M+NH4-H2O]+', '[M+H-C13H12O9]+', 
                             '[M+H-C3H8NO6P]+', '[M+H-CH4O]+', '[M+H-C2H6O]+', '[M+H-CH4]+', '[M+H-C11H12N2O3]+', '[M-H-C10H20]-', 
                             '[M+CH3OH+H]+', '[M+H-C8H10O]+', '[M-C3H7O2]-', '[M-H-C3H5NO2]-', '[M+H-C5H14NO4P]+', '[M+H-C24H44O-H2O]+', 
                             '[M+H-C12H20O9]+', '[M+H-C6H10O5]+', '[M+H-C6H10O4]+', '[M-CH3]-', '[M-C2H3O]-', '[M-C3H8O+H]+', '[M+H-C5H12N2]+',
                              '[M+H-C4H6]+', '[M-H-C6H9O5SO3H]-', '[M+H-C2H5N]+', '[M-H-CO2]-', '[M+H-2CH4]+', '[M+H-C5H9NO4]+', '[M+H-CH3NH2]+', 
                              '[M-H-CO2-2HF]-', '[M+H-3CH4]+',  '[M+K-2H]-', '[M-H2+H]+', '[M-SO3+H]+', '[M-H-NH3]-']

metadata = metadata[~metadata['ADDUCT'].isin(remove_neutralloss_adduct)]

remove_multimer_adduct = ['[2M-3H2O]+', '[2M-2H2O+NH4]+',  '[2M-2H2O]+', '[2M-3H2O+NH4]+', '[2M-H2O+NH4]+', '[2M-H2O]+', '[2M]+', 
                          '[2M-3H2O+Na]+', '[2M-H2O+Na]+', '[2M+K]+', '[2M+NH4]+', '[2M+ACN+H]+', '[2M+Hac-H]-', '[2M+K-2H]-', 
                          '[4M-H]-', '[2M-2H+Na]-', '[2M-2H+K]-','[M+3ACN+2H]2+', '[M+2Na]2+', '[M-2H]2-', '[M-3H]3-', '[M+H+K]2+', 
                          '[M+H+NH4]2+', '[M+2ACN+2H]2+', '[M+2NH4]2+']

metadata = metadata[~metadata['ADDUCT'].isin(remove_multimer_adduct)]

unknown_adduct = ['[M-e]-', '[M+Li]+', '[M-H+Li]+', '[M+Br]-', '[M+TFA-H]-', '[M+MeOH-H]-',  '[M+2K-H]+', '[M+IsoProp+H]+', '[M+2Na-H]+', 
                  '[M+DMSO+H]+', '[M+ACN+H]+', '[M+ACN+Na]+', '[2M-2H+3Na]+'] 

metadata = metadata[~metadata['ADDUCT'].isin(unknown_adduct)]

additional_adduct = ['[M-2H]-', '[M-H]+', '[M+H]-', '[M+2H]+', '[M+CH3]+', '[M+OH]+', '[M+H+O]+', '[M-H+H2O]-', '[M+H+H2O]+', '[M+Na-2H]-', 
                     '[M-H+Na]+', '[M+H+Na]+', '[M-H+CH3OH]-', '[M+K]-', '[M+C2H3N+NH4]+', '[M+OAc]-', '[M+FA+H]+']

metadata = metadata[~metadata['ADDUCT'].isin(additional_adduct)]

metadata['ADDUCT'] = metadata['ADDUCT'].replace('[M+Hac-H]-', '[M+CH3COO]-')

metadata['ppm_error'] = ''

# ref: https://github.com/matchms/matchms/blob/1f904e0d469aef35dbba8b7b2d7b52886f3f75cc/matchms/data/known_adducts_table.csv#L56
adduct_mass_adjustments = {
    '[M]+': 0,
    '[M+K]+': 38.963158,
    '[M+CH3COO]-': 59.013851,
    '[M+H+HCOOH]+': 46.00548,
    '[M]-': 0,
    '[M+FA-H]-': (46.005477 - 1.007276),
    '[M+Cl]-' : 34.969402,
    '[M+NH4]+': 18.033823,
    '[M+Na]+' : 22.989218,
    '[M-H]-': -1.007276,
    '[M+H]+': 1.007276,
}

def ppm_error(expt, theo):
    return abs((expt - theo) / expt) * 1e6

for i in range(len(metadata)):
    adduct = metadata['ADDUCT'].iloc[i]
    if adduct in adduct_mass_adjustments:
        theo = metadata['PARENT_MASS'].iloc[i] + adduct_mass_adjustments[adduct]
        expt = metadata['PRECURSOR_MZ'].iloc[i]
        metadata['ppm_error'].iloc[i] = ppm_error(expt, theo)     

metadata = metadata[metadata['ppm_error'] <= 10]

adduct_m_minus = metadata[metadata['ADDUCT'] == '[M]-']

def calculate_properties(smile):
    mol = Chem.MolFromSmiles(smile)
    if mol:
        formula = rdMolDescriptors.CalcMolFormula(mol)
        exact_mass = rdMolDescriptors.CalcExactMolWt(mol)
        charge = Chem.GetFormalCharge(mol)
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.InchiToInchiKey(inchi)
        return formula, charge, exact_mass, inchi, inchikey
    return None, None, None, None, None

for index, row in adduct_m_minus.iterrows():
    if row['SMILES'] == 'CS(=O)(=O)CCCC(=NOS(=O)(=O)[O-])SC1OC(CO)C(O)C(O)C1O':
        smile = 'CS(=O)(CCCC(SC1OC(C(C(C1O)O)O)CO)=NOS(=O)(O)=O)=O'
    elif row['SMILES'] == 'O=S(=O)([O-])ON=C(Cc1ccccc1)SC1OC(CO)C(O)C(O)C1O':
        smile = 'O=S(O)(ON=C(SC1OC(C(C(C1O)O)O)CO)Cc2ccccc2)=O'
    elif row['SMILES'] == 'O=S(=O)([O-])ON=C(Cc1ccc(O)cc1)SC1OC(CO)C(O)C(O)C1O':
        smile = 'O=S(O)(ON=C(SC1OC(C(C(C1O)O)O)CO)Cc2ccc(O)cc2)=O'
    elif row['SMILES'] == 'CC(C)CCCC(C)C1CCC2C3CC=C4CC(OS(=O)(=O)[O-])CCC4(C)C3CCC12C':
        smile = 'CC(CCCC(C1CCC2C3CC=C4CC(CCC4(C3CCC12C)C)OS(=O)(O)=O)C)C'
    elif row['SMILES'] == 'CC(C)(O)C(=O)[O-]':
        smile = 'CC(O)(C(O)=O)C'
    else:
        raise ValueError(f"Unknown SMILES: {row['SMILES']}")

    formula, charge, exact_mass, inchi, inchikey = calculate_properties(smile)
    metadata.at[index, 'SMILES'] = smile
    metadata.at[index, 'INCHIKEY'] = inchikey
    metadata.at[index, 'FORMULA'] = formula
    metadata.at[index, 'CHARGE'] = charge
    metadata.at[index, 'PARENT_MASS'] = exact_mass
    metadata.at[index, 'INCHI'] = inchi
    metadata.at[index, 'ADDUCT'] = '[M-H]-'

adduct_m_minus = metadata[metadata['ADDUCT'] == '[M]-']

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

with open(output_mgf_dir, 'w') as file:
    for sublist in info:
        for item in sublist:
            file.write("%s\n" % item)
        file.write("\n") 

metadata.to_csv(output_csv_dir, index=False)