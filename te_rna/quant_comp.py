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

parser = argparse.ArgumentParser(description="Compare TE expression using BLAST + Salmon TPMs per condition")

parser.add_argument("--blast", required=True,
                    help="BLAST tabular file (outfmt 6)")
parser.add_argument("--tpm_instances", required=True,
                    help="Salmon merged TPM file for TE instances")
parser.add_argument("--tpm_representatives", required=True,
                    help="Salmon merged TPM file for representative TEs")
parser.add_argument("--condition_map", required=True,
                    help="TSV file mapping sample replicates to biological condition (columns: sample, condition)")
parser.add_argument("--min_cov", type=float, default=0.8,
                    help="Minimum coverage required (default = 0.8 = 80%)")
parser.add_argument("--out", default="TE_expression_comparison.tsv",
                    help="Output comparison table")

args = parser.parse_args()

# Load condition map

print(f"Using condition map from file: {args.condition_map}")
cond_df = pd.read_csv(args.condition_map, sep="\t")
condition_map = dict(zip(cond_df["sample"], cond_df["condition"]))

# Mean TPMs per condition

def mean_tpm_per_condition(tpm_df, condition_map):
    cond_to_samples = defaultdict(list)
    for sample, cond in condition_map.items():
        if sample in tpm_df.columns:
            cond_to_samples[cond].append(sample)

    cond_mean = pd.DataFrame(index=tpm_df.index)

    for cond, samples in cond_to_samples.items():
        cond_mean[cond] = tpm_df[samples].mean(axis=1)
    
    return cond_mean

# Load data from BLAST

blast_cols = [
    "instance",
    "representative",
    "pident",
    "aln_len",
    "qlen",
    "slen",
    "evalue",
    "bitscore"
]

blast = pd.read_csv(
    args.blast,
    sep="\t",
    header=None,
    names=blast_cols
)

# Filter BLAST by coverage (double check from file already filtered)

print(f"Filtering BLAST hits with coverage >= {args.min_cov*100:.0f}%")

blast["cov_q"] = blast["aln_len"] / blast["qlen"]
blast_filt = blast[blast["cov_q"] >= args.min_cov]

# Map instances to representatives

# mapping = defaultdict(list)
mapping = defaultdict(set)

for _, row in blast_filt.iterrows():
#   mapping[row["representative"]].append(row["instance"])
    mapping[row["representative"]].add(row["instance"])

# Load TPM tables

print("Loading TPM tables...")

tpm_inst = pd.read_csv(args.tpm_instances, sep="\t", index_col=0)
tpm_repr = pd.read_csv(args.tpm_representatives, sep="\t", index_col=0)

tpm_inst = mean_tpm_per_condition(tpm_inst, condition_map)
tpm_repr = mean_tpm_per_condition(tpm_repr, condition_map)

# Consistency checks:
# Representatives that exist structurally but show no expression
structural_only = set(mapping.keys()) - set(tpm_repr.index)

print(f"{len(structural_only)} representatives TEs detected by BLAST, "
      f"but with no detectable expression")

# Expressed TEs without BLAST support (QC issue)
expression_only = set(tpm_repr.index) - set(mapping.keys())
if expression_only:
    print(f"Warning: {len(expression_only)} expressed TEs "
          f"not supported by BLAST mapping")

# Sum instance TPMs per TE, per condition

inst_sum = defaultdict(dict)

for te, inst_list in mapping.items():
    for cond in tpm_repr.columns:
        total = 0.0
        for inst in inst_list:
            if inst in tpm_inst.index:
                total += tpm_inst.loc[inst, cond]
        inst_sum[te][cond] = total
        
# Build comparison table

rows = []

for te in tpm_repr.index:
    for cond in tpm_repr.columns:
        rows.append({
            "TE": te,
            "condition": cond,
            "TPM_representative": tpm_repr.loc[te, cond],
            "TPM_instances_sum": inst_sum.get(te, {}).get(cond, 0),
            "n_instances": len(mapping.get(te, []))
    })

out_df = pd.DataFrame(rows)
out_df.to_csv(args.out, sep="\t", index=False)

print(f"Comparison written to {args.out}")
