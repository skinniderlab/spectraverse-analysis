import pandas as pd
import os, sys

print("Step1-2: MGF to CSV conversion")

mgf_dir = sys.argv[1]
csv_dir = sys.argv[2] 

os.makedirs(mgf_dir, exist_ok=True)
os.makedirs(csv_dir, exist_ok=True)

files = os.listdir(mgf_dir)

for file in files:
    filename = os.path.join(mgf_dir, file)
    print(filename)

    rows = []
    row = {}

    with open(filename, 'r') as file:
        for line in file:
            line = line.rstrip('\n')

            if line == 'END IONS':
                rows.append(row)
                row = {}
            elif '=' in line and not line.split('=')[0].strip().isdigit():
                column, value = line.split('=', 1)
                row[column] = value

    base_name = os.path.basename(filename)
    new_name = base_name.replace('.mgf', '.csv')
    df = pd.DataFrame(rows)
    assert not df['SOURCE'].isna().any(), 'check if any row for a column is NaN'
    df.to_csv(os.path.join(csv_dir, new_name), index=False)




