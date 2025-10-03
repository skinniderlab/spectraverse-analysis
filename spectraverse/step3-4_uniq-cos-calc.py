import numpy as np
import os, glob
import pandas as pd
import spectral_entropy
import sys
from matchms.importing import load_from_mgf
from matchms.filtering import normalize_intensities
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm
from multiprocessing import Pool

print("Step3-4: Calculating pairwise cosine scores for duplicate spectra and saving in numpy format")

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

create_dir_and_delete_files(numpy_dir)

dedup = metadata.copy()

dedup['INCHI14'] = dedup['INCHIKEY'].str[:14]
dedup = dedup.assign(unique_identifier = dedup[['INCHI14',
                                                'IONMODE']].astype(str).\
                     agg('_'.join, axis=1))

dup = dedup[dedup.duplicated(subset=['unique_identifier'])]
dup_identifiers = dup['unique_identifier'].unique()
dup_df = dedup[dedup['unique_identifier'].isin(dup_identifiers)]

# define function to perform dot_product
def pairwise_dot(spectra_arr, out):
    pairwise_dot_nb(spectra_arr, out)
    # Guard against numerical instability.
    return np.clip(out, 0, 1, out)

def pairwise_dot_nb(spectra_arr, out):
    for i in range(len(spectra_arr)):
        for j in range(i + 1, len(spectra_arr)):
            out[i, j] = out[j, i] = dot(spectra_arr[i], spectra_arr[j])
    np.fill_diagonal(out, 1)

def dot(spectrum_query, spectrum_library):
    score = spectral_entropy.similarity(
                spectrum_query, 
                spectrum_library, 
                method='dot_product',
                need_clean_spectra=False, 
                ms2_ppm = 5
                )
    return score 

# calculate dot_prodcuts for files with same identifier and save it in numpy format
def process_identifier(identifier):
    mgf_filepath = mgf_dir + '/{}.mgf'.format(identifier)
    numpy_filepath = numpy_dir + '/{}.npy'.format(identifier)
    
    # normalize intensities
    spectra = []
    for spectrum in load_from_mgf(mgf_filepath):
        spectrum = normalize_intensities(spectrum)
        spectra.append(spectrum)
    
    spectra_list = []
    for i in range(len(spectra)):
        mz = spectra[i].mz.reshape(1, -1)
        if mz.shape[1] == 1:
            mz = np.c_[mz, 0]
        intensities = spectra[i].intensities.reshape(1, -1)
        if intensities.shape[1] == 1:
            intensities = np.c_[intensities, 0]
        mz_intensity = np.concatenate((mz, intensities), 0)
        mz_intensity = np.transpose(mz_intensity)
        spectra_list.append(mz_intensity)
    
    dot_products = pairwise_dot(spectra_list, 
                                np.zeros((len(spectra_list), 
                                          len(spectra_list)), np.float32))
    np.save(numpy_filepath, dot_products)

with Pool() as pool:
    list(tqdm(pool.imap(process_identifier, dup_identifiers), total=len(dup_identifiers)))