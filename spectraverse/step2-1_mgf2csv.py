import pandas as pd
import os, sys

print("Step2-1: MGF to CSV conversion")

csv_dir = sys.argv[1]
mgf_dir = sys.argv[2]

rows = [] 
row = {}

with open(mgf_dir, 'r') as file:
    for line in file:
        line = line.rstrip('\n')

        if line == 'END IONS':
            rows.append(row)
            row = {}
        elif '=' in line and not line.split('=')[0].strip().isdigit():
            column, value = line.split('=', 1)
            row[column] = value

df = pd.DataFrame(rows)
df.to_csv(csv_dir, index=False)




