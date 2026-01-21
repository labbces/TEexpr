"""
Compare TE expression using representative TEs vs genomic instances.

Steps:
1. Read BLAST results (instances vs representatives)
2. Filter alignments by minimum coverage (default: 80%)
3. Map instances to representative TEs
4. Sum TPMs of instances per representative
5. Compare with TPMs from representative-based quantification
"""

import argparse
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser(description="Compare TE expression using BLAST + Salmon TPMs")

parser.add_argument("--blast", required=True, help="BLAST tabular file (outfmt 6)")
parser.add_argument("--tpm_instances", required=True, help="Salmon merged TPM file for TE instances")
parser.add_argument("--tpm_representatives", required=True, help="Salmon merged TPM file for representative TEs")
parser.add_argument("--min_cov", type=float, default=0.8, help="Minimum coverage required (default = 0.8 = 80%)")
parser.add_argument("--out", default="TE_expression_comparison.tsv", help="Output comparison table")

args = parser.parse_args()

# -----------------------------
# LOAD DATA
# -----------------------------

print("Loading BLAST results...")
blast = pd.read_csv(
    args.blast,
    sep="\t",
    header=None,
    names=["instance", "representative", "pident",
           "aln_len", "qlen", "slen", "evalue", "bitscore"]
)

print("Loading Salmon TPMs...")
tpm_inst = pd.read_csv(args.tpm_instances, sep="\t", index_col=0)
tpm_repr = pd.read_csv(args.tpm_representatives, sep="\t", index_col=0)

# -----------------------------
# FILTER BLAST BY COVERAGE
# -----------------------------

print(f"Filtering BLAST hits with coverage >= {args.min_cov*100:.0f}%")

blast["cov_q"] = blast["aln_len"] / blast["qlen"]
blast_filt = blast[blast["cov_q"] >= args.min_cov]

# -----------------------------
# MAP INSTANCES TO REPRESENTATIVES
# -----------------------------

mapping = defaultdict(list)

for _, row in blast_filt.iterrows():
    mapping[row["representative"]].append(row["instance"])

# -----------------------------
# SUM INSTANCE TPMs PER TE
# -----------------------------

inst_sum = {}

for te, inst_list in mapping.items():
    total = 0.0
    for inst in inst_list:
        if inst in tpm_inst.index:
            total += tpm_inst.loc[inst].sum()
    inst_sum[te] = total

# -----------------------------
# BUILD COMPARISON TABLE
# -----------------------------

rows = []

for te in tpm_repr.index:
    rows.append({
        "TE": te,
        "TPM_representative": tpm_repr.loc[te].sum(),
        "TPM_instances_sum": inst_sum.get(te, 0),
        "n_instances": len(mapping.get(te, []))
    })

out_df = pd.DataFrame(rows)
out_df.to_csv(args.out, sep="\t", index=False)

print(f"Comparison written to {args.out}")
