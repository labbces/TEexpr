from collections import Counter

from Bio import Entrez
import pandas as pd
from io import StringIO
import time
import argparse

def parse_args():
	parser = argparse.ArgumentParser(
		description="Retrieve SRRs from BioProjects and compare with a local list"
	)

	parser.add_argument(
		"-b",
		"--bioprojects",
		required=True,
		help="CSV file containing BioProject IDs"
	)

	parser.add_argument(
		"-s",
		"--species",
		required=True,
		help="Species name (e.g., 'Homo sapiens')"
	)

	parser.add_argument(
		"-e",
		"--email",
		required=True,
		help="Email for Entrez (NCBI requirement)"
	)

	parser.add_argument(
		"-k",
		"--api_key",
		help="NCBI API Key for increased request limits",
		required=True
	)

	parser.add_argument(
		"-o",
		"--output",
		default="sra_report",
		help="Output file prefix"
	)

	return parser.parse_args()


# =======================================

def load_bioprojects(csv_path):
	"""
	Expects a CSV with a column containing BioProject IDs.
	Accepts:
	- column named 'bioproject'
	- otherwise uses the first column
	"""
	df = pd.read_csv(csv_path)

	if "bioproject" in df.columns:
		return df["bioproject"].dropna().astype(str).unique().tolist()
	else:
		return df.iloc[:, 0].dropna().astype(str).unique().tolist()


def fetch_sra_runinfo(bioproject):
	print(f"[INFO] Fetching {bioproject}...")

	handle = Entrez.esearch(
		db="sra",
		term=bioproject,
		retmax=100000
	)
	record = Entrez.read(handle)
	ids = record["IdList"]

	if not ids:
		print(f"[WARNING] No SRA records found for {bioproject}")
		return pd.DataFrame()

	print(f"[INFO] Retrieved {len(ids)} records for {bioproject}")

	handle = Entrez.efetch(
		db="sra",
		id=",".join(ids),
		rettype="runinfo",
		retmode="text"
	)

	data = handle.read().decode('utf-8')  # Decodifica bytes para string
	df = pd.read_csv(StringIO(data))

	df["BioProject"] = bioproject

	time.sleep(0.12)  # avoid NCBI rate limits
	return df


def load_srr_list(csv_path):
	df = pd.read_csv(csv_path)
	if "srr" in df.columns:
		return set(df["srr"].dropna().astype(str))
	else:
		return set(df.iloc[:, 0].dropna().astype(str))


def analyze_bioprojects(bioproject_list, species):
	all_dfs = []

	for bp in bioproject_list:
		df = fetch_sra_runinfo(bp)
		if not df.empty:
			all_dfs.append(df)

	if not all_dfs:
		raise ValueError("No data retrieved from NCBI")

	df_total = pd.concat(all_dfs, ignore_index=True)

	# More robust species filtering (handles strain/isolate variations)
	df_species = df_total[
		df_total["ScientificName"].str.contains(species, na=False)
	]

	return df_total, df_species


def generate_report(df_species, user_srrs, output_prefix, species, bioproject_list):
	srrs_expected = set(df_species["Run"])
	srrs_user = user_srrs

	missing = srrs_expected - srrs_user
	extra = srrs_user - srrs_expected

	print("\n===== SUMMARY =====")
	print(f"Total SRRs (species): {len(srrs_expected)}")
	print(f"In your list: {len(srrs_user)}")
	print(f"Missing: {len(missing)}")
	print(f"Extra: {len(extra)}")

	# Save filtered metadata
	df_species.to_csv(f"{output_prefix}_filtered.csv", index=False)

	# Save text report
	with open(f"{output_prefix}.txt", "w") as f:
		f.write("===== SRA REPORT =====\n\n")
		f.write(f"Species: {species}\n")
		f.write(f"BioProjects: {', '.join(bioproject_list)}\n\n")

		f.write(f"Total expected SRRs: {len(srrs_expected)}\n")
		f.write(f"SRRs in your list: {len(srrs_user)}\n\n")

		f.write("=== MISSING SRRs ===\n")
		for s in sorted(missing):
			f.write(s + "\n")

		f.write("\n=== EXTRA SRRs ===\n")
		for s in sorted(extra):
			f.write(s + "\n")
		
		f.write("\n=== DUPLICATED SRRs IN YOUR LIST ===\n")

	print(f"\n[OK] Reports saved with prefix: {output_prefix}")


# =======================================

def main():
	args = parse_args()

	Entrez.email = args.email
	Entrez.api_key = args.api_key

	bioproject_list = load_bioprojects(args.bioprojects)

	user_srrs = load_srr_list(args.bioprojects)
	
	df_total, df_species = analyze_bioprojects(
		bioproject_list,
		args.species
	)

	generate_report(
		df_species,
		user_srrs,
		args.output,
		args.species,
		bioproject_list
	)


if __name__ == "__main__":
	main()