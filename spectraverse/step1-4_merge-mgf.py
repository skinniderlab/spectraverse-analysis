import os, sys
import glob
import pandas as pd
import numpy as np

print("Step1-4: Merge MGF files")

mgf_dir = sys.argv[1]
output_file = sys.argv[2]

os.makedirs(os.path.dirname(output_file), exist_ok=True)

mgf_files = [f for dirpath, dirnames, filenames in os.walk(mgf_dir) for f in glob.glob(os.path.join(dirpath, '*.mgf'))]
mgf_files = sorted(mgf_files)
mgf_file_names = [os.path.basename(f) for f in mgf_files]

with open(output_file, 'w') as outfile:
    for file_name in mgf_files:
        print(file_name)
        with open(file_name, 'r') as infile:
            lines = infile.readlines()

            while lines and not lines[-1].strip():
                lines.pop()

            outfile.write(''.join(lines) + '\n\n')