import os, sys
import glob
import pandas as pd
import numpy as np

print("Step1-3: Merge CSV files")

csv_dir = sys.argv[1]
output_file = sys.argv[2]

os.makedirs(os.path.dirname(output_file), exist_ok=True)

csv_files = [f for dirpath, dirnames, filenames in os.walk(csv_dir) for f in glob.glob(os.path.join(dirpath, '*.csv'))]
csv_files = sorted(csv_files)
csv_file_names = [os.path.basename(f) for f in csv_files]

csv_dfs = []
columns = ['INDEX', 'TITLE', 'SOURCE' ,'SMILES', 'INCHI','INCHIKEY', 'INCHIAUX', 'IONMODE', 'PRECURSOR_MZ', 'PEPMASS', 'CHARGE', 'NUM_PEAKS', 'ADDUCT', 'NUM_PEAKS', 'MS_LEVEL', 'COMPOUND_NAME', 'SPECTRUMID']

for csv_file in csv_files:
    print(csv_file)
    df = pd.read_csv(csv_file)
    for col in columns:
        if col not in df.columns:
            df[col] = ''

    df_temp = df[columns]
    csv_dfs.append(df_temp)

csv_merge = pd.concat(csv_dfs, ignore_index=True)

# Correct adducts   
csv_merge.loc[csv_merge['ADDUCT'] == '[M+FA]-', 'ADDUCT'] = '[M+FA-H]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+NH4]-', 'ADDUCT'] = '[M+NH4]+'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+H2CO2-H]1-', 'ADDUCT'] = '[M+HCO2]-'
csv_merge.loc[csv_merge['ADDUCT'] == '(M+H)+', 'ADDUCT'] = '[M+H]+'
csv_merge.loc[csv_merge['ADDUCT'] == '[M-HCl+H]+', 'ADDUCT'] = '[M-Cl]+'
csv_merge.loc[csv_merge['ADDUCT'] == 'M+', 'ADDUCT'] = '[M]+'
csv_merge.loc[csv_merge['ADDUCT'] == '[M-H+OH]-', 'ADDUCT'] = '[M+O]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+HCOO]-', 'ADDUCT'] = '[M+HCO2]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[Cat]+', 'ADDUCT'] = '[M]+'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+C2H3N+Na-2H]-', 'ADDUCT'] = '[M+C2H3N+Na-2H]-'
csv_merge.loc[csv_merge['ADDUCT'] == 'M + Formate', 'ADDUCT'] = '[M+C2F3O2]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+]', 'ADDUCT'] = '[M]+'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+Formate]', 'ADDUCT'] = '[M+C2F3O2]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+acetate]', 'ADDUCT'] = '[M+CH3COO]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+formate]', 'ADDUCT'] = '[M+C2F3O2]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M-]', 'ADDUCT'] = '[M]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+Cl-]', 'ADDUCT'] = '[M+Cl]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+FA-]', 'ADDUCT'] = '[M+FA-H]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M-H1]', 'ADDUCT'] = '[M-H]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[(M+CH3COOH)-H]-', 'ADDUCT'] = '[M+CH3COOH-H]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M-1]-', 'ADDUCT'] = '[M-H]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+C]', 'ADDUCT'] = '[M+Cl]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+HAc-H]-', 'ADDUCT'] = '[M+CH3COO]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+H-SO3]+', 'ADDUCT'] = '[M-SO3+H]+'
csv_merge.loc[csv_merge['ADDUCT'] == '[M-NH4+2H]+', 'ADDUCT'] = '[M+2H-NH4]+'
csv_merge.loc[csv_merge['ADDUCT'] == 'Deprotonated molecule', 'ADDUCT'] = '[M-H]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+CH3COOH-H]-', 'ADDUCT'] = '[M+CH3COO]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M-H]1-', 'ADDUCT'] = '[M-H]-'
csv_merge.loc[csv_merge['ADDUCT'] == '[M+ACN+NH4]+', 'ADDUCT'] = '[M+C2H3N+NH4]+'
csv_merge.loc[csv_merge['ADDUCT'] == '[M-e]', 'ADDUCT'] = '[M-e]-'

matchms_adduct_list = ['[M+Na-2H]', '[M+H+Na]+', 'M+Na-2H', '[M-HCl-H]-', '[M-2HCl+H]', '[M-HCl+Na]+', '[M-HCl+K]+', '[M+NH3]', 
               '[M-Ac-H-]', '[?]', 'carotenoid', 'carotenoids', '[M-H+C2H2O]', '[Unknown]', '[Cat]', 'MSMS', '[Unk]', '[--]']

csv_merge['ADDUCT'] = csv_merge['ADDUCT'].apply(lambda x: None if x in matchms_adduct_list else x)

csv_merge.to_csv(output_file, index=False)




