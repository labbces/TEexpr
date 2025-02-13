#!/usr/bin/env python3

import os
import pandas as pd

def process_files(species_list, output_dir):
    # Garantir que o diretório de saída existe
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    for species in species_list:
        # Definir caminhos dos arquivos
        earl_grey_file = f"{species}_earl_grey.fasta"
        deep_te_file = f"{species}_deep_te.fasta"
        te_sorter_file = f"{species}_te_sorter.cls.tsv"
        output_file = os.path.join(output_dir, f"{species}_comparativo.csv")
        
        # Processar Earl Grey
        with open(earl_grey_file, 'r') as f:
            earl_grey_headers = [line.strip()[1:] if line.startswith('>') else line.strip() for line in f if line.startswith('>')]
        
        # Processar DeepTE
        with open(deep_te_file, 'r') as f:
            deep_te_headers = {}
            for line in f:
                if line.startswith('>'):
                    parts = line.strip().split("__")
                    if len(parts) > 1:
                        deep_te_headers[parts[0][1:]] = parts[1]
                    else:
                        deep_te_headers[parts[0][1:]] = "Unknown"
        
        # Processar TEsorter
        te_sorter_df = pd.read_csv(te_sorter_file, sep='\t')
        columns = te_sorter_df.columns

        if len(columns) < 4:
            raise ValueError(f"O arquivo {te_sorter_file} não possui as colunas esperadas.")

        te_sorter_classes = {}
        for _, row in te_sorter_df.iterrows():
            id_part = row[columns[0]].split('#')[0]  # Extrai o ID antes do '#'
            combined_classes = "_".join(map(str, row[columns[1:4]].tolist()))
            te_sorter_classes[id_part] = combined_classes

        # Criar tabela final
        with open(output_file, 'w') as out:
            out.write("Earl_Grey,DeepTE,TEsorter\n")
            for header in earl_grey_headers:
                deep_te_value = deep_te_headers.get(header, "Not_Found")
                te_sorter_value = te_sorter_classes.get(header, "Not_Found")
                out.write(f"{header},{deep_te_value},{te_sorter_value}\n")
    
    print("Processamento concluído! Arquivos salvos em:", output_dir)

# Exemplo de uso:
species_list = ["species1", "species2", "species3"]  # Substitua pelas suas espécies
output_dir = "output"
process_files(species_list, output_dir)

