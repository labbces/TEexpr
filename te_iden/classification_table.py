#!/usr/bin/env python3

import os
import pandas as pd

def process_files(species_list, output_dir):
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    for species in species_list:
        # Define file paths
        earl_grey_file = f"{species}_earl_grey.fasta"
        deep_te_file = f"{species}_deep_te.fasta"
        te_sorter_file = f"{species}_te_sorter.cls.tsv"
        output_file = os.path.join(output_dir, f"{species}_comparativo.csv")
        
        # Process Earl Grey
        with open(earl_grey_file, 'r') as f:
            earl_grey_headers = [line.strip()[1:] if line.startswith('>') else line.strip() for line in f if line.startswith('>')]
        
        # Process DeepTE
        with open(deep_te_file, 'r') as f:
            deep_te_headers = {}
            for line in f:
                if line.startswith('>'):
                    parts = line.strip().split("__")
                    if len(parts) > 1:
                        deep_te_headers[parts[0][1:]] = parts[1]
                    else:
                        deep_te_headers[parts[0][1:]] = "Unknown"
        
        # Process TEsorter
        te_sorter_df = pd.read_csv(te_sorter_file, sep='\t')
        columns = te_sorter_df.columns

        if len(columns) < 4:
            raise ValueError(f"The file {te_sorter_file} does not have the expected columns.")

        te_sorter_classes = {}
        for _, row in te_sorter_df.iterrows():
            id_part = row[columns[0]].split('#')[0]  # Extract the ID before '#'
            combined_classes = "_".join(map(str, row[columns[1:4]].tolist()))
            te_sorter_classes[id_part] = combined_classes

        # Create final table
        with open(output_file, 'w') as out:
            out.write("Earl_Grey,DeepTE,TEsorter\n")
            for header in earl_grey_headers:
                deep_te_value = deep_te_headers.get(header, "Not_Found")
                te_sorter_value = te_sorter_classes.get(header, "Not_Found")
                out.write(f"{header},{deep_te_value},{te_sorter_value}\n")
    
    print("Processing completed! Files saved in:", output_dir)

# Usage example:
species_list = ["species1", "species2", "species3"]  # Replace with your species
output_dir = "output"
process_files(species_list, output_dir)

