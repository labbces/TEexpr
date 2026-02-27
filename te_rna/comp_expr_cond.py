"""
Compare TE expression using representative TEs vs genomic instances.
Steps:
1. Read RepeatMasker output (instances vs representatives)
2. Filter alignments by minimum identity (default: 80%)
3. Map instances to representative TEs
4. Average TPMs of instances per representative, per condition
5. Compare with TPMs mean from representative-based quantification, per condition
"""

import argparse
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser(description="Compare TE expression using RepeatMasker + Salmon TPMs")
parser.add_argument("--rm", required=True, help="RepeatMasker .out file")
parser.add_argument("--tpm_inst", required=True, help="Salmon merged TPM file for TE instances")
parser.add_argument("--tpm_repr", required=True, help="Salmon merged TPM file for representative TEs")
parser.add_argument("--cond", required=True, help="TSV with columns: sample, condition")
parser.add_argument("--min_identity", type=float, default=80.0, help="Minimum identity required (default = 80.0%%)")
parser.add_argument("--out", default="TE_expression_comparison.tsv", help="Output comparison table")
args = parser.parse_args()

# -----------------------------
# LOAD DATA
# -----------------------------
print("Loading RepeatMasker output...")
rm_cols = ["SW_score", "perc_div", "perc_del", "perc_ins",
           "query_seq", "q_start", "q_end", "q_left", "strand",
           "matching_repeat", "repeat_class",
           "r_start", "r_end", "r_left", "ID"]

rm = pd.read_csv(
    args.rm,
    sep=r"\s+",
    header=None,
    skiprows=2,
    names=rm_cols,
    on_bad_lines="skip"
)

# Compute identity and filter
rm["identity"] = 100 - rm["perc_div"].astype(float)
rm_filt = rm[rm["identity"] >= args.min_identity]
print(f"  {len(rm)} total hits → {len(rm_filt)} kept after identity >= {args.min_identity}%")

print("Loading Salmon TPMs...")
tpm_inst = pd.read_csv(args.tpm_inst, sep="\t", index_col=0)
tpm_repr = pd.read_csv(args.tpm_repr, sep="\t", index_col=0)

print("Loading conditions...")
cond_df = pd.read_csv(args.cond, sep="\t")
sample_to_cond = dict(zip(cond_df["sample"], cond_df["condition"]))

# -----------------------------
# MAP INSTANCES TO REPRESENTATIVES
# -----------------------------
# query_seq = genomic instance; matching_repeat = representative TE
mapping = defaultdict(set)
for _, row in rm_filt.iterrows():
    mapping[row["matching_repeat"]].add(row["query_seq"])

# -----------------------------
# GROUP SAMPLES BY CONDITION
# -----------------------------
cond_to_samples = defaultdict(list)
for sample, cond in sample_to_cond.items():
    if sample in tpm_inst.columns:
        cond_to_samples[cond].append(sample)

missing = [s for s in sample_to_cond if s not in tpm_inst.columns]
if missing:
    print(f"Warning: these samples are in the conditions file but not in tpm_inst: {missing}")

# -----------------------------
# BUILD COMPARISON TABLE
# (one row per TE per condition)
# -----------------------------
rows = []
for cond, samples in cond_to_samples.items():
    inst_cond = tpm_inst[samples]
    repr_cond = tpm_repr[[s for s in samples if s in tpm_repr.columns]]

    for te in tpm_repr.index:
        tpm_rep_val = repr_cond.loc[te].mean() if te in repr_cond.index else 0.0

        tpm_inst_val = 0.0
        for inst in mapping.get(te, []):
            if inst in inst_cond.index:
                tpm_inst_val += inst_cond.loc[inst].mean()

        rows.append({
            "TE":                 te,
            "condition":          cond,
            "TPM_representative_mean": tpm_rep_val,
            "TPM_instances_mean":  tpm_inst_val,
            "n_instances":        len(mapping.get(te, []))
        })

out_df = pd.DataFrame(rows)
out_df.to_csv(args.out, sep="\t", index=False)
print(f"Comparison written to {args.out}")
