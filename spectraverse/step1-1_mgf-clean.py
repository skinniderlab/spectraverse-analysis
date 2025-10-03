import os, sys, re
import glob
import pandas as pd
import numpy as np

print("Step1-1: MGF cleaning and metadata standardization")

input_file_dir = sys.argv[1]
output_file_dir = sys.argv[2]

GNPS_FILE_NAME = 'test_gnps.mgf'
MONA_FILE_NAME = 'test_mona.mgf'
MASSBANK_FILE_NAME = 'test_massbank.mgf'

os.makedirs(output_file_dir, exist_ok=True)

mgf_files = [f for dirpath, dirnames, filenames in os.walk(input_file_dir) for f in glob.glob(os.path.join(dirpath, '*.mgf'))]
input_file_name = [os.path.basename(f) for f in mgf_files]
output_file_name = input_file_name 

for i in range(len(input_file_name)):
    source_text = 'SOURCE={}\n'.format(input_file_name[i])

    count = 0
    index = 0
    inside_block = False
    block_lines = []
    with open(input_file_dir+input_file_name[i], 'r') as source, open(output_file_dir+output_file_name[i], 'w') as destination:
        for line in source:
            if "BEGIN IONS" in line:
                inside_block = True
                block_lines.append(line)
                block_lines.append(source_text)
                block_lines.append("INDEX={}\n".format(index))
                if input_file_name[i] != GNPS_FILE_NAME or input_file_name[i] != MONA_FILE_NAME or input_file_name[i] != MASSBANK_FILE_NAME:     
                    block_lines.append("MS_LEVEL=2\n")
                precursor_mz = 0
                pepmass = 0
                num_peaks = 0
                index = index + 1
            elif inside_block and "END IONS" in line:
                inside_block = False                      
                if not any("NUM_PEAKS" in block_line for block_line in block_lines) or (count!=num_peaks):
                    block_lines.insert(1, f"NUM_PEAKS={count}\n")
                if not any("PEPMASS=" in block_line for block_line in block_lines) and precursor_mz!=0:
                    block_lines.insert(1, f"PEPMASS={precursor_mz}\n") 
                if not any("PRECURSOR_MZ=" in block_line for block_line in block_lines) and pepmass!=0:
                    block_lines.insert(1, f"PRECURSOR_MZ={pepmass}\n")   
                if not any("MS_LEVEL=" in block_line for block_line in block_lines):       
                    if input_file_name[i] != GNPS_FILE_NAME or input_file_name[i] != MONA_FILE_NAME or input_file_name[i] != MASSBANK_FILE_NAME: 
                        block_lines.insert(1, "MS_LEVEL=2\n")
                block_lines.append(line)
                for block_line in block_lines:
                    destination.write(block_line)
                count = 0
                block_lines = []
            elif inside_block and line[0].isdigit():
                parts = line.split()
                if float(parts[1]) != 0.0 and float(parts[0]) > 10.0 and float(parts[0]) < 1000.0:
                    count += 1
                    line = ' '.join(parts) + '\n'
                    block_lines.append(line)
            elif line.startswith("NUM_PEAKS"):
                match = re.search(r'NUM_PEAKS=(.*)', line)
                if match and match.group(1).strip() != "":
                    num_peaks = int(match.group(1).strip())         
            elif line.startswith("CHARGE"):    
                match = re.search(r'CHARGE=(.*)', line)
                if match and match.group(1).strip() != "":
                    block_lines.append(line)
            elif line.startswith("charge"):    
                match = re.search(r'charge=(.*)', line)
                if match and match.group(1).strip() != "":
                    line = line.replace("charge", "CHARGE")
                    block_lines.append(line)
            elif line.startswith("MSLEVEL"):    
                match = re.search(r'MSLEVEL=(.*)', line)
                if match and match.group(1).strip() != "":
                    line = line.replace("MSLEVEL", "MS_LEVEL")
                    block_lines.append(line)      
            elif line.startswith("mslevel"):    
                match = re.search(r'mslevel=(.*)', line)
                if match and match.group(1).strip() != "":
                    line = line.replace("mslevel", "MS_LEVEL")
                    block_lines.append(line)            
            elif line.startswith("PRECURSOR_MZ"):
                match = re.search(r'PRECURSOR_MZ=(.*)', line)
                if match and match.group(1).strip() != "":
                    precursor_mz = float(match.group(1).strip())
                    block_lines.append(line)  
            elif line.startswith("PEPMASS"):
                match = re.search(r'PEPMASS=(.*)', line)
                if match and match.group(1).strip() != "":
                    pepmass = float(match.group(1).strip())
                    block_lines.append(line)
            elif line.startswith("IONMODE"):
                match = re.search(r'IONMODE=(.*)\n', line)
                if match and match.group(1).strip() != "":
                    ionmode = match.group(1).strip()
                    if ionmode in ['p', 'Positive', 'POS', 'e', 'positive']:
                        block_lines.append("IONMODE=positive\n")
                    elif ionmode in ['n', 'Negative', 'NEG', 'negative']:
                        block_lines.append("IONMODE=negative\n")                  
            # GNPS Specific preprocessing
            elif inside_block and "NAME=" in line and "FILENAME=" not in line:
                if input_file_name[i] == GNPS_FILE_NAME:  
                    adduct = line.rsplit(' ', 1)[-1]
                    adduct = adduct.replace("\n", "")
                    if '[' in adduct:
                        block_lines.append('ADDUCT=' + adduct +'\n')
                    else:    
                        block_lines.append('ADDUCT=[' + adduct +']\n')    
            # In-house libraires specific preprocessing         
            elif line.startswith("METABOLITE_IDENTIFICATION"):    
                if input_file_name[i] != GNPS_FILE_NAME or input_file_name[i] != MONA_FILE_NAME or input_file_name[i] != MASSBANK_FILE_NAME:  
                    match = re.search(r'METABOLITE_IDENTIFICATION=(.*)', line)
                    if match and match.group(1).strip() != "":
                        line = line.replace("METABOLITE_IDENTIFICATION", "COMPOUND_NAME")
                        block_lines.append(line) 
            elif line.startswith("MODIFICATIONS"):    
                if input_file_name[i] != GNPS_FILE_NAME or input_file_name[i] != MONA_FILE_NAME or input_file_name[i] != MASSBANK_FILE_NAME:  
                    match = re.search(r'MODIFICATIONS=(.*)', line)
                    if match and match.group(1).strip() != "":
                        line = line.replace("MODIFICATIONS", "ADDUCT")
                        block_lines.append(line)  

            elif inside_block:
                block_lines.append(line)
            else:
                destination.write(line)