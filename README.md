# <img align="left" width="80" height="80" style="vertical-align: bottom;" src="assets/logo.png" alt="Spectraverse icon"> Comprehensive curation and harmonization of small molecule MS/MS libraries in Spectraverse

This repository contains code used to preprocess and harmonize MS/MS spectra and their associated metadata in the compilation of Spectraverse.

## Environment Setup

### Prerequisites

- Python 3.10 or higher
- pip or conda package manager

### Option 1: Using Conda (Recommended)

1. **Install Conda**: Download and install [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

2. **Create Environment**:
   ```bash
   conda env create -f environment.yml
   conda activate spectraverse
   ```

3. **Verify Installation**:
   ```bash
   python --version
   pip list
   ```

### Option 2: Using Virtual Environment

1. **Create Virtual Environment**:
   ```bash
   python3.10 -m venv spectraverse
   source spectraverse/bin/activate 
   ```

2. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/GuptaVishu2002/spectraverse.git
   cd spectraverse
   ```

2. **Set up Environment** (choose one of the options above)


3. **Install Package in Development Mode**:
   ```bash
   pip install -e .
   ```

4. **Install Other Package**:
     ```bash
   # Code used to perform matchms preprocessing   
   git clone https://github.com/matchms/matchms
   # To use our version of matchms libraries with modified code use the zip file at
   # matchms/matchms.zip and replace it original matchms libraries
   # this version fixes the function _get_neutral_mass to ignore spectra annotated
   # with [M]+ or [M]− adducts that had a reported charge of zero.
     

   # Code used to calculate cosine score
   git clone https://github.com/YuanyueLi/SpectralEntropy.git
   cd SpectralEntropy
   mv spectral_entropy spectraverse  # Move spectral_entropy folder inside the spectraverse directory
   ```
   


## Usage

### Quick Start

1. **Prepare your data**: Place raw data files (MGF format) in the `data/` directory.
2. **Configure processing**: Edit `config/config_step1.json`, `config/config_step2.json` and `config/config_step3.json` (see JSON files for documentation).
3. **Run the master script**: Execute the preprocessing pipeline by running `run_steps.py`.
4. **Check results**: Find processed data in `data/`

### Example Usage

The main script `run_steps.py` orchestrates the entire preprocessing pipeline. It provides a unified interface to run all preprocessing steps in sequence. 

```bash

# Step 1: Run the initial preprocessing.
python run_steps.py --config config/config_step1.json

# Step 2: Run matchms using filters specified in matchms/library_cleaning_1.yaml.

# Step 3:Run the second preprocessing step.
python run_steps.py --config config/config_step2.json

# Step 4: Run matchms again using filters specified in matchms/library_cleaning_2.yaml.

# Step 5: Run the final preprocessing step.
python run_steps.py --config config/config_step3.json

```
Note: It is crucial to execute the steps in this order, as each stage depends on the outputs of the previous step.


### Configuration

The config file contains the following keys to specify python script and data paths.

- `script_base_path`: Specify python script path
- `data_base_path`: Specify data file path
- `pipeline`: Specify the order of python scripts to run
- `steps`: Specify initial, intermidiate and final data file names

### Pipeline Steps

The main script executes the following steps in order:

1. `config_step1.json`
- **step1-1_mgf-clean.py**: MGF cleaning and metadata standardization
- **step1-2_mgf2csv.py**: MGF to CSV conversion
- **step1-3_merge-csv.py**: Merge CSV files
- **step1-4_merge-mgf.py**: Merge MGF files
- **step1-5_csv2mgf.py**: CSV to MGF conversion (including some manual standardization of metadata fields)

At this point, matchms is run to repair the metadata associated with each spectrum.

2. `config_step2.json`
- **step2-1_mgf2csv.py**: MGF to CSV conversion
- **step2-2_cleaning.py**: Cleaning metadata and MGF files by removing spectra with identical fragment intensities
- **step2-3_standardization.py**: Standardize SMILES
- **step2-4_removal.py**: Removal of unwanted spectra based on various criteria
- **step2-5_lowres-check.py**: Removal of low-resolution spectra
- **step2-6_prec-check.py**: Removal of spectra based on precursor and fragment mass check
- **step2-7_highfragmass-check.py**: Removal of spectra with all fragment mass > PRECURSOR_MZ
- **step2-8_csv2mgf.py**: CSV to MGF conversion (accomplished with metadata modification)

At this point, matchms is run to again repair metadata and remove incoherent annotations.

3. `config_step3.json`
- **step3-1_mgf2csv.py**: MGF to CSV conversion
- **step3-2_adduct-rem.py**: Removal and correction of unwanted adducts
- **step3-3_uniq-comb.py**: Identification of candidate duplicate spectra
- **step3-4_uniq-cos-calc.py**: Calculating pairwise cosine scores for duplicate spectra and saving in numpy format
- **step3-5_uniq-select.py**: Removal of duplicate spectra based on pairwise cosine scores
- **step3-6_noise-rem.py**: Removal of low-intensity fragment ions and additional structurally uninformative spectra
- **step3-7_metadata.py**: Standardization of instrument and collision energy fields, and finalizing the metadata
