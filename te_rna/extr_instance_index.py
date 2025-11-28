#!/usr/bin/env python3
"""
Script otimizado para extrair cópias (instâncias) genômicas de TEs.

Requisitos:
- Python 3.8+
- biopython
- pyranges

Instalar dependências:
	pip install biopython pyranges
"""
import argparse
import sys
import re
import pandas as pd
import pyranges as pr
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parse_repeatmasker_out(rm_path):
	#Parseia arquivo .out do RepeatMasker e retorna DataFrame.
	rows = []
	with open(rm_path) as fh:
		for raw in fh:
			line = raw.rstrip('\n')
			if not line:
				continue

			# Ignora cabeçalhos do RepeatMasker
			if line.lstrip().startswith('SW') or line.lstrip().startswith('score'):
				continue

			fields = line.split()

			# Linhas de hits começam com número (score)
			if re.match(r'^\d+$', fields[0]):
				try:
					# Extrai informações principais do hit
					contig = fields[4]
					start = int(fields[5]) - 1  # RM é base-1, Python é base-0
					end = int(fields[6])
					strand = '+'

					# Determina a fita (strand) do TE
					if fields[8] == 'C' or fields[8] == 'c':
						strand = '-'
						repeat_field_idx = 9
					else:
						if fields[8] in ('+', '-'):
							strand = fields[8]
							repeat_field_idx = 9
						else:
							repeat_field_idx = 8
					name = fields[repeat_field_idx]

					# Adiciona ao resultado
					rows.append({
						'Chromosome': contig,
						'Start': start,
						'End': end,
						'Strand': strand,
						'TE_name': name
					})
				except Exception:
					# Se não conseguir parsear, ignora a linha
					continue

	# Retorna DataFrame para facilitar manipulação posterior
	return pd.DataFrame(rows)

def parse_bed(bed_path):
	#Parseia arquivo BED e retorna DataFrame.
	rows = []
	with open(bed_path) as fh:
		for raw in fh:
			line = raw.strip()
			if not line or line.startswith('#'):
				continue
			parts = line.split('\t')
			if len(parts) < 3:
				continue

			# Extrai informações do BED
			contig = parts[0]
			start = int(parts[1])
			end = int(parts[2])
			name = parts[3] if len(parts) > 3 else 'TE'
			strand = parts[5] if len(parts) > 5 else '+'
			rows.append({
				'Chromosome': contig,
				'Start': start,
				'End': end,
				'Strand': strand,
				'TE_name': name
			})

	# Retorna DataFrame para facilitar manipulação posterior
	return pd.DataFrame(rows)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Extrai cópias genômicas de TEs (RepeatMasker .out ou BED) e gera FASTA para Salmon usando PyRanges e Biopython.')
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument('--repeatmasker', '-r', help='RepeatMasker .out file')
	group.add_argument('--bed', '-b', help='BED file with TE positions')
	parser.add_argument('--genome', '-g', required=True, help='Genome FASTA')
	parser.add_argument('--out_fasta', '-o', required=True, help='Output FASTA with TE genomic instances')
	parser.add_argument('--min_length', type=int, default=20, help='Minimum length to keep an instance')
	parser.add_argument('--filter', '-f', help='Regex: only include TE names matching this pattern')
	parser.add_argument('--no_simplify_names', action='store_true', help='Do not simplify TE names (keep original characters)')
	args = parser.parse_args()

	# Carregar hits do arquivo de entrada
	if args.repeatmasker:
		print('Parsing RepeatMasker .out...')
		df = parse_repeatmasker_out(args.repeatmasker)
	else:
		print('Parsing BED...')
		df = parse_bed(args.bed)

	# Filtrar por regex (nome do TE)
	if args.filter:
		df = df[df['TE_name'].str.contains(args.filter, regex=True)]
		print(f"Filtered hits; {len(df)} remaining after applying pattern {args.filter}")

	# Filtrar por tamanho mínimo de instância
	df = df[(df['End'] - df['Start']) >= args.min_length]

	# Simplificar nomes para evitar problemas em identificadores
	if not args.no_simplify_names:
		df['TE_name'] = df['TE_name'].str.replace(r'[^A-Za-z0-9_.-]', '_', regex=True)

	# Criar objeto PyRanges para manipulação eficiente dos intervalos
	gr = pr.PyRanges(df)

	# Carregar genoma em dicionário de cromossomos
	print('Carregando genoma...')
	genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))

	# Inicializa listas para registros e metadados
	records = []
	tsv_lines = []
	map_lines = []

	# Contador de cópias por família de TE
	te_counts = defaultdict(int)
	# Itera sobre cada intervalo (instância de TE)
	for row in gr.df.itertuples(index=False):
		te_counts[row.TE_name] += 1
		copy_num = te_counts[row.TE_name]

		# Gera identificador único para cada cópia
		copy_id = f"{row.TE_name}_copy{copy_num:04d}|{row.Chromosome}:{row.Start+1}-{row.End}|{row.Strand}"

		try:
			# Extrai sequência do genoma
			seq_obj = genome[row.Chromosome].seq[row.Start:row.End]
		except KeyError:
			# Se cromossomo não existe, ignora
			print(f"Contig {row.Chromosome} não encontrado no genoma. Pulando {copy_id}", file=sys.stderr)
			continue

		# Aplica reverse complement se necessário
		if row.Strand == '-':
			seq_obj = seq_obj.reverse_complement()

		seq = str(seq_obj)
		seqlen = len(seq)
		# Ignora instâncias curtas
		if seqlen < args.min_length:
			continue
		# Cria registro FASTA e metadados
		records.append(SeqRecord(seq_obj, id=copy_id, description=f"family={row.TE_name} location={row.Chromosome}:{row.Start+1}-{row.End}({row.Strand})"))
		tsv_lines.append(f"{copy_id}\t{row.TE_name}\t{row.Chromosome}\t{row.Start+1}\t{row.End}\t{row.Strand}\t{seqlen}")
		map_lines.append(f"{copy_id}\t{row.TE_name}")

	# Escrever arquivos de saída (FASTA, TSV, MAP)
	print(f"Escrevendo {len(records)} sequências em {args.out_fasta}...")
	SeqIO.write(records, args.out_fasta, "fasta")

	with open(args.out_fasta + '.tsv', 'w') as f:
		f.write('\t'.join(['copy_id', 'TE_name', 'contig', 'start', 'end', 'strand', 'length']) + '\n')
		for line in tsv_lines:
			f.write(line + '\n')

	with open(args.out_fasta + '.transcript_to_gene_map.tsv', 'w') as f:
		for line in map_lines:
			f.write(line + '\n')
	print('Arquivos gerados com sucesso.')

	# Comandos sugeridos para uso com Salmon
	salmon_index_cmd = f"salmon index -t {args.out_fasta} -i salmon_index_TEgenomic -k 31 --transcript-to-gene-map {args.out_fasta}.transcript_to_gene_map.tsv"
	salmon_quant_cmd = "salmon quant -i salmon_index_TEgenomic -l A -1 <reads_R1.fastq.gz> -2 <reads_R2.fastq.gz> -o <sample_TEgenomic_quant>"

	print('\nSuggested Salmon commands:')
	print(salmon_index_cmd)
	print(salmon_quant_cmd)
