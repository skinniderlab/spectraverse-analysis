import numpy as np
import os, glob, re
import pandas as pd
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdMolDescriptors
from tqdm import tqdm
from multiprocessing import Pool

print("Step2-3: Standardize SMILES")

metadata_csv_dir = sys.argv[1]
input_mgf_dir = sys.argv[2]

benchmark_dir = os.path.dirname(metadata_csv_dir)

temp_mgf_dir = benchmark_dir + '/test-temp.mgf'
temp_csv_dir = benchmark_dir + '/test-temp.csv'

output_csv_dir = sys.argv[3]
output_mgf_dir = sys.argv[4]


removed_index = []

metadata = pd.read_csv(metadata_csv_dir)

deuterated = metadata[metadata['SMILES'].fillna('').str.contains('[2H]', regex=False)].index
metadata = metadata[~metadata['SMILES'].fillna('').str.contains('[2H]', regex=False)]

tritiated = metadata[metadata['SMILES'].fillna('').str.contains('[3H]', regex=False)].index
metadata = metadata[~metadata['SMILES'].fillna('').str.contains('[3H]', regex=False)]

c13_smiles = metadata[metadata['SMILES'].fillna('').str.contains('[13C]', regex=False)].index
metadata = metadata[~metadata['SMILES'].fillna('').str.contains('[13C]', regex=False)]

precursor_mz_1000 = metadata[metadata['PRECURSOR_MZ'] > 1000]
metadata = metadata[metadata['PRECURSOR_MZ'] <= 1000]

metadata_columns = list(metadata.columns)

metadata_no_smiles = metadata[metadata['SMILES'].isna()]
no_smiles_dir = '/Genomics/argo/users/vg8892/git/msms-triangulation/data/spectra/benchmark/no_smiles.csv'
metadata_no_smiles.to_csv(no_smiles_dir, index=False)

metadata = metadata[metadata['SMILES'].notna()]
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

with open(temp_mgf_dir, 'w') as file:
    for sublist in info:
        for item in sublist:
            file.write("%s\n" % item)
        file.write("\n") 

metadata.to_csv(temp_csv_dir, index=False)

"""
Step 3: standardize SMILES and add inchikeys
"""
uncharger = rdMolStandardize.Uncharger()
te = rdMolStandardize.TautomerEnumerator()

# define a function to check for aromatic sulfoxides and correct bond/charges
def check_aromatic_sulfoxides(parent_clean_mol):
    S_pos = []
    O_pos = []
    S_charge = False
    O_charge = False
    single_bond = False
    for atom in parent_clean_mol.GetAtoms():
        if atom.GetSymbol() == 'S':
            charge = atom.GetFormalCharge()
            if charge == 1:
                S_charge = True
                S_pos.append(atom.GetIdx())
        if atom.GetSymbol() == 'O':
            charge = atom.GetFormalCharge()
            degree = atom.GetDegree()
            if charge == -1 and degree == 1:
                O_charge = True    
                O_pos.append(atom.GetIdx())
    bonds = []
    SO_comb = []
    SO_comb_idx = []
    for i in range(len(S_pos)):
        for j in range(len(O_pos)):
            bonds.append(parent_clean_mol.GetBondBetweenAtoms(S_pos[i], O_pos[j]))
            SO_comb.append((S_pos[i], O_pos[j]))
    none_indices = [i for i, x in enumerate(bonds) if x is None]
    bonds = [x for x in bonds if x is not None]
    SO_comb = [x for i, x in enumerate(SO_comb) if i not in none_indices]        
    for i in range(len(bonds)):
        if bonds[i].GetBondType() == Chem.rdchem.BondType.SINGLE:
            single_bond = True
            SO_comb_idx.append(i)
    
    if S_charge and O_charge and single_bond:
        for i in range(len(SO_comb_idx)):
            atom1 = parent_clean_mol.GetAtomWithIdx(SO_comb[SO_comb_idx[i]][0])
            atom1.SetFormalCharge(0)
            atom2 = parent_clean_mol.GetAtomWithIdx(SO_comb[SO_comb_idx[i]][1])
            atom2.SetFormalCharge(0)
            bond = parent_clean_mol.GetBondBetweenAtoms(SO_comb[SO_comb_idx[i]][0], SO_comb[SO_comb_idx[i]][1])
            bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
    
    return parent_clean_mol

def check_P_O_charge(parent_clean_mol):
    P_pos = []
    O_pos = []
    P_charge = False
    O_charge = False
    single_bond = False
    for atom in parent_clean_mol.GetAtoms():
        if atom.GetSymbol() == 'P':
            charge = atom.GetFormalCharge()
            if charge == 1:
                P_charge = True
                P_pos.append(atom.GetIdx())
        if atom.GetSymbol() == 'O':
            charge = atom.GetFormalCharge()
            degree = atom.GetDegree()
            if charge == -1 and degree == 1:
                O_charge = True    
                O_pos.append(atom.GetIdx())
    bonds = []
    PO_comb = []
    PO_comb_idx = []
    for i in range(len(P_pos)):
        for j in range(len(O_pos)):
            bonds.append(parent_clean_mol.GetBondBetweenAtoms(P_pos[i], O_pos[j]))
            PO_comb.append((P_pos[i], O_pos[j]))
    none_indices = [i for i, x in enumerate(bonds) if x is None]
    bonds = [x for x in bonds if x is not None]
    PO_comb = [x for i, x in enumerate(PO_comb) if i not in none_indices]
    for i in range(len(bonds)):
        if bonds[i].GetBondType() == Chem.rdchem.BondType.SINGLE:
            single_bond = True
            PO_comb_idx.append(i)

    if P_charge and O_charge and single_bond:
        for i in range(len(PO_comb_idx)):
            atom1 = parent_clean_mol.GetAtomWithIdx(PO_comb[PO_comb_idx[i]][0])
            atom1.SetFormalCharge(0)
            atom2 = parent_clean_mol.GetAtomWithIdx(PO_comb[PO_comb_idx[i]][1])
            atom2.SetFormalCharge(0)
            bond = parent_clean_mol.GetBondBetweenAtoms(PO_comb[PO_comb_idx[i]][0], PO_comb[PO_comb_idx[i]][1])
            bond.SetBondType(Chem.rdchem.BondType.DOUBLE)

    return parent_clean_mol

"""
rdkit contributed code to neutralize charged molecules;
obtained from:
    https://www.rdkit.org/docs/Cookbook.html
    http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg02669.html
"""
def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None
def NeutraliseCharges(mol, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    return mol

def preprocess_mol(smiles, stereochem=False, standardize_tautomers=True):
    """
    SMILES preprocessing and standardization pipeline:
    1. Load molecule in RDKit
    2. Optionally remove stereochemistry
    3. Sanitize and remove hydrogens
    4. Run 'Cleanup' (disconnect metal atoms, reionize)
    5. Select parent molecule if multiple fragments
    6. Neutralize charges, using two different approaches
    7. Optionally, standardize tautomers
    8. Return the cleaned molecule, canonical SMILES, and inchikey.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # raise ValueError("invalid SMILES: " + str(smiles))
        print("invalid SMILES: " + str(smiles))
        return None, None, None, None, None, None, None, None 
    
    # optionally, remove stereochemistry
    if not stereochem:
        Chem.RemoveStereochemistry(mol)
    
    # sanitize molecule and remove hydrogens
    try:
    # sanitize molecule and remove hydrogens
        Chem.SanitizeMol(mol)
    except:
        print("Unable to kekulize molecule: " + str(smiles))
        return None, None, None, None, None, None, None, None 
    mol = Chem.RemoveHs(mol)
    
    # from https://bitsilla.com/blog/2021/06/standardizing-a-molecule-using-rdkit/
    # removeHs, disconnect metal atoms, normalize the molecule, reionize the molecule
    try:
        clean_mol = rdMolStandardize.Cleanup(mol) 
    except:
        print("Unable to clean molecule: " + str(smiles))
        return None, None, None, None, None, None, None, None    
    # if many fragments, get the "parent" (the actual mol we are interested in) 
    try:
        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)
    except:
        print("Unable to fragment parent molecule: " + str(smiles))
        return None, None, None, None, None, None, None, None
    
    parent_clean_mol = check_aromatic_sulfoxides(parent_clean_mol)
    parent_clean_mol = check_P_O_charge(parent_clean_mol)
    
    charge_before = sum(atom.GetFormalCharge() for atom in parent_clean_mol.GetAtoms())
    
    # two different approaches to neutralizing charges
    uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)
    # manual double-check
    uncharged_parent_clean_mol2 = NeutraliseCharges(uncharged_parent_clean_mol)
    # run uncharger again
    uncharged_parent_clean_mol3 = uncharger.uncharge(uncharged_parent_clean_mol2)
    
    charge_after = sum(atom.GetFormalCharge() for atom in uncharged_parent_clean_mol3.GetAtoms())
    
    # tautomer enumerator
    if standardize_tautomers:
        try:
            Chem.SanitizeMol(uncharged_parent_clean_mol3)
            mol = te.Canonicalize(uncharged_parent_clean_mol3)
        except:
            print("Unable to canonicalize tautomers: " + str(smiles))
            mol = smiles = inchikey = None   
    else:
        mol = uncharged_parent_clean_mol3
        try:
            Chem.SanitizeMol(mol)
        except:
            print("Unable to sanitize molecule: " + str(smiles))
            mol = smiles = inchikey = None    
    
    """
    Not used: 
        Chem.GetSymmSSSR(mol)  # forces RDKit to find rings and perceive aromaticity
        Chem.SanitizeMol(mol, Chem.SANITIZE_SETAROMATICITY | Chem.SANITIZE_KEKULIZE)
    """
    
    # get SMILES and inchikey
    if mol is not None:
        try:
            smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=stereochem)

            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.RemoveHs(mol)
            smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=stereochem)

            inchikey = Chem.inchi.MolToInchiKey(mol)
        except:
            print("Unable to convert molecule to SMILES or InChIKey: " + str(smiles))
            smiles = inchikey = None
    else:
        smiles = inchikey = None
    
    if charge_before != charge_after:
        duplicate = 'duplicate'
        # tautomer enumerator
        if standardize_tautomers:
            try:
                Chem.SanitizeMol(parent_clean_mol)
                mol_charge = te.Canonicalize(parent_clean_mol)
            except:
                print("Unable to canonicalize tautomers: " + str(smiles))
                mol_charge = smiles_charge = inchikey_charge = None    
        else:
            mol_charge = parent_clean_mol
            try:
                Chem.SanitizeMol(mol_charge)
            except:
                print("Unable to sanitize molecule: " + str(smiles))
                mol_charge = smiles_charge = inchikey_charge = None  
        
        if mol_charge is not None:
            try:
                smiles_charge = Chem.MolToSmiles(mol_charge, canonical=True, isomericSmiles=stereochem)

                mol_charge = Chem.MolFromSmiles(smiles_charge)
                mol_charge = Chem.RemoveHs(mol_charge)
                smiles_charge = Chem.MolToSmiles(mol_charge, canonical=True, isomericSmiles=stereochem)

                inchikey_charge = Chem.inchi.MolToInchiKey(mol_charge)  
            except:
                print("Unable to convert molecule to SMILES or InChIKey: " + str(smiles))
                smiles_charge = inchikey_charge = None    
    else:
        duplicate = 'unique'
        mol_charge = smiles_charge = inchikey_charge = None
    
    return mol, smiles, inchikey, mol_charge, smiles_charge, inchikey_charge, duplicate

def process_smiles(smiles):
    try:
        smiles = smiles.split(' |')[0]
        mol, canonical_sm, inchikey, mol_charge, canonical_sm_charge, inchikey_charge, duplicate = preprocess_mol(smiles)
        return mol, canonical_sm, inchikey, mol_charge, canonical_sm_charge, inchikey_charge, duplicate
    except ValueError as e:
        return None, None, None, None, None, None, None

metadata = pd.read_csv(temp_csv_dir)

uniq_smiles = metadata.drop_duplicates('SMILES')

with Pool() as p:
    results = p.map(process_smiles, uniq_smiles['SMILES'].values)

mols, canonical_smiles, inchikeys, charge_cols, charge_canonical_smiles, charge_inchikey, duplicate = zip(*results)
# add to data frame
uniq_smiles = uniq_smiles.assign(CANONICAL_SMILES=canonical_smiles, CANONICAL_INCHIKEY=inchikeys,
                                 CANONICAL_SMILES_CHARGE=charge_canonical_smiles, CANONICAL_INCHIKEY_CHARGE=charge_inchikey,
                                 DUPLICATE=duplicate)
uniq_smiles = uniq_smiles[['SMILES', 'CANONICAL_SMILES', 'CANONICAL_INCHIKEY', 'CANONICAL_SMILES_CHARGE', 'CANONICAL_INCHIKEY_CHARGE', 'DUPLICATE']]
metadata = metadata.merge(uniq_smiles, how='left', on='SMILES')

metadata1 = metadata[metadata_columns + ['CANONICAL_SMILES', 'CANONICAL_INCHIKEY', 'DUPLICATE']]
metadata2 = metadata[metadata_columns + ['CANONICAL_SMILES_CHARGE', 'CANONICAL_INCHIKEY_CHARGE', 'DUPLICATE']]

metadata2 = metadata2.rename(columns={'CANONICAL_SMILES_CHARGE': 'CANONICAL_SMILES', 'CANONICAL_INCHIKEY_CHARGE': 'CANONICAL_INCHIKEY'})

metadata1 = metadata1[~metadata1['CANONICAL_SMILES'].isna()]
metadata2 = metadata2[~metadata2['CANONICAL_SMILES'].isna()]

metadata = pd.concat([metadata1, metadata2])
metadata = metadata.sort_index()
metadata_index = metadata.index

def calculate_exact_mass(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        exact_mass = Descriptors.ExactMolWt(mol)
    except:
        exact_mass = None    

    return exact_mass

metadata['EXACT_MASS_ORIGINAL'] = metadata['SMILES'].apply(calculate_exact_mass)
metadata['EXACT_MASS_CANONICAL'] = metadata['CANONICAL_SMILES'].apply(calculate_exact_mass)

for i in range(metadata.shape[0]):
    smiles = metadata['CANONICAL_SMILES'].iloc[i]
    mol = Chem.MolFromSmiles(smiles)
    charge = Chem.GetFormalCharge(mol)
    if charge != 0:
        parent_mass = metadata['PARENT_MASS'].iloc[i]
        exact_mass = metadata['EXACT_MASS_CANONICAL'].iloc[i]
        if metadata['ADDUCT'].iloc[i] == '[M+H]+':
            if abs(parent_mass - exact_mass) <= 0.1 and charge == 1:
                metadata['ADDUCT'].iloc[i] = '[M]+'
        if metadata['ADDUCT'].iloc[i] == '[M-H]-':
            if abs(parent_mass - exact_mass) <= 0.1 and charge == -1:
                metadata['ADDUCT'].iloc[i] = '[M]-'

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
    '[M]': 0,
    '[M+K]': 38.963158,
    '[M+CH3COO]': 59.013851,
    '[M+H+HCOOH]': 46.00548,
    '[M+FA-H]': (46.005477 - 1.007276),
    '[M+Cl]' : 34.969402,
    '[M+NH4]': 18.033823,
    '[M+Na]' : 22.989218,
    '[M-H]': -1.007276,
    '[M+H]': 1.007276,
}

def normalize_adduct(adduct):
    if pd.isna(adduct):
        return adduct    
    return adduct.rstrip('+-')

for i in range(metadata.shape[0]):
    orig_adduct = metadata['ORIG_ADDUCT'].iloc[i]
    adduct = metadata['ADDUCT'].iloc[i]
    if normalize_adduct(orig_adduct) != normalize_adduct(adduct) and orig_adduct in adduct_mass_adjustments and adduct in adduct_mass_adjustments:
        parent_mass = metadata['EXACT_MASS_CANONICAL'].iloc[i]
        precusor_mz = metadata['PRECURSOR_MZ'].iloc[i]
    
        mass_diff = precusor_mz - parent_mass

        orig_adduct_value = adduct_mass_adjustments.get(orig_adduct)
        adduct_value = adduct_mass_adjustments.get(adduct)

        orig_adduct_diff = abs(mass_diff - orig_adduct_value)
        adduct_diff = abs(mass_diff - adduct_value)

        if orig_adduct_diff < adduct_diff:
            metadata['ADDUCT'].iloc[i] = orig_adduct

info = extract_info(temp_mgf_dir, metadata_index)

with open(output_mgf_dir, 'w') as file:
    for sublist in info:
        for item in sublist:
            file.write("%s\n" % item)
        file.write("\n") 

metadata.to_csv(output_csv_dir, index=False)

