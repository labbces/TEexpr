#!/usr/bin/env python3
"""
Compare TE expression using representative TEs vs genomic instances.
Steps:
1. Read metadata TSV from extract_te_copies.py (copy_id to TE_name mapping)
2. Map instances to representative TEs
3. Average TPMs of instances per representative, per condition
4. Compare with TPMs mean from representative-based quantification, per condition

Usage example:
    python comp_expr_cond.py \\
        --metadata   EDTA_instances.fasta.tsv \\
        --tpm_instances      EDTA_tpm_instances.tsv \\
        --tpm_representatives EDTA_tpm_representatives.tsv \\
        --conditions conditions.tsv \\
        --out        EDTA_quant_comp.tsv
"""
import argparse
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser(description="Compare TE expression using RepeatMasker + Salmon TPMs")
parser.add_argument("--metadata",            required=True,  help="Metadata TSV from extract_te_copies.py (out_fasta.tsv)")
parser.add_argument("--tpm_instances",       required=True,  help="Salmon merged TPM file for TE instances")
parser.add_argument("--tpm_representatives", required=True,  help="Salmon merged TPM file for representative TEs")
parser.add_argument("--conditions",          required=True,  help="TSV with columns: sample, condition")
parser.add_argument("--out",  default="TE_expression_comparison.tsv", help="Output comparison table")
args = parser.parse_args()

# -----------------------------
# LOAD DATA
# -----------------------------
print("Loading metadata from extract_te_copies.py...")
meta = pd.read_csv(args.metadata, sep="\t")
# Expected columns: copy_id, TE_name, contig, start, end, strand, length
print(f"  {len(meta)} instances loaded.")

print("Loading Salmon TPMs...")
tpm_inst = pd.read_csv(args.tpm_instances, sep="\t", index_col=0)
tpm_repr = pd.read_csv(args.tpm_representatives, sep="\t", index_col=0)

print("Loading conditions...")
cond_df = pd.read_csv(args.conditions, sep="\t")
sample_to_cond = dict(zip(cond_df["sample"], cond_df["condition"]))

# -----------------------------
# MAP INSTANCES TO REPRESENTATIVES
# using copy_id (matches tpm_inst index) → TE_name (matches tpm_repr index)
# -----------------------------
mapping = defaultdict(set)
for _, row in meta.iterrows():
    mapping[row["TE_name"]].add(row["copy_id"])

# Sanity check: how many copy_ids are found in tpm_inst
all_copy_ids = {cid for cids in mapping.values() for cid in cids}
found = all_copy_ids & set(tpm_inst.index)
print(f"  {len(found)} of {len(all_copy_ids)} instances found in tpm_instances.")

# -----------------------------
# GROUP SAMPLES BY CONDITION
# -----------------------------
cond_to_samples = defaultdict(list)
for sample, cond in sample_to_cond.items():
    if sample in tpm_inst.columns:
        cond_to_samples[cond].append(sample)

missing = [s for s in sample_to_cond if s not in tpm_inst.columns]
if missing:
    print(f"Warning: these samples are in the conditions file but not in tpm_instances: {missing}")

# -----------------------------
# BUILD COMPARISON TABLE
# (one row per TE per condition)
# -----------------------------
rows = []
for cond, samples in cond_to_samples.items():
    inst_cond = tpm_inst[samples]
    repr_cond = tpm_repr[[s for s in samples if s in tpm_repr.columns]]

    for te in tpm_repr.index:
        te_id = te.split("#")[0]  # strip #CLASS → "TE_00000778"
        tpm_rep_val = repr_cond.loc[te].mean() if te in repr_cond.index else 0.0

        tpm_inst_val = 0.0
        for inst in mapping.get(te_id, []):
            if inst in inst_cond.index:
                tpm_inst_val += inst_cond.loc[inst].mean()

        rows.append({
            "TE":                      te,
            "condition":               cond,
            "TPM_representative_mean": tpm_rep_val,
            "TPM_instances_mean":      tpm_inst_val,
            "n_instances":             len(mapping.get(te_id, []))
        })

out_df = pd.DataFrame(rows)
out_df.to_csv(args.out, sep="\t", index=False)
print(f"Comparison written to {args.out}")
