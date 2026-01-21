#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 10
#$ -l h=neotera

# Load environments
module load blast/2.16.0+
module load Python/3.12.2
source /Storage/data2/andreza.cunha/scripts/te_instance_env/bin/activate

# -----------------------------
# PATHS
# -----------------------------

BASE="/Storage/data1/andreza.cunha/rna_seq/salmon_test/comparative_test"
SCRIPT_COMPARE="/Storage/data2/andreza.cunha/scripts/quant_comp.py"

# -----------------------------
# EARL GREY FILES
# -----------------------------

EG_REPR="/Storage/data1/andreza.cunha/rna_seq/salmon_test/index_test/distilledTE.flTE.iter655.fa"
EG_INST="/Storage/data1/andreza.cunha/rna_seq/salmon_test/index_test/EG_maize_instances.fa"
EG_DB="$BASE/eg_blast_db"
EG_BLAST="$BASE/blast_EG_comp.tsv"

EG_TPM_INST="$BASE/maize_EG_instances_out/merged_TPM.tsv"
EG_TPM_REPR="$BASE/maize_EG_distilled_out/merged_TPM.tsv"
EG_OUT_COMPARE="$BASE/EG_inst_vs_repr_TPM.tsv"

# -----------------------------
# EDTA FILES
# -----------------------------

EDTA_REPR="/Storage/data1/andreza.cunha/rna_seq/salmon_test/index_test/distilledTE.flTE.iter3274.fa"
EDTA_INST="/Storage/data1/andreza.cunha/rna_seq/salmon_test/index_test/EDTA_maize_instances.fa"
EDTA_DB="$BASE/edta_blast_db"
EDTA_BLAST="$BASE/blast_EDTA_inst_vs_repr.tsv"

EDTA_TPM_INST="$BASE/maize_EDTA_instances_out/merged_TPM.tsv"
EDTA_TPM_REPR="$BASE/maize_EDTA_distilled_out/merged_TPM.tsv"
EDTA_OUT_COMPARE="$BASE/EDTA_inst_vs_repr_TPM.tsv"

# -----------------------------
# EARL GREY — BLAST
# -----------------------------

echo ">>> BLAST — Earl Grey"
echo "Using $NSLOTS CPUs"

makeblastdb -in "$EG_REPR" -dbtype nucl -out "$EG_DB"

blastn \
  -query "$EG_INST" \
  -db "$EG_DB" \
  -num_threads $NSLOTS \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  -out "$EG_BLAST"

# -----------------------------
# EARL GREY — TPM COMPARISON
# -----------------------------

echo ">>> Comparing Salmon TPMs — Earl Grey"

python3 "$SCRIPT_COMPARE" \
  --blast "$EG_BLAST" \
  --tpm_instances "$EG_TPM_INST" \
  --tpm_representatives "$EG_TPM_REPR" \
  --min_cov 0.8 \
  --out "$EG_OUT_COMPARE"

# -----------------------------
# EDTA — BLAST
# -----------------------------

echo ">>> BLAST — EDTA"
echo "Using $NSLOTS CPUs"

makeblastdb -in "$EDTA_REPR" -dbtype nucl -out "$EDTA_DB"

blastn \
  -query "$EDTA_INST" \
  -db "$EDTA_DB" \
  -num_threads $NSLOTS \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  -out "$EDTA_BLAST"

# -----------------------------
# EDTA — TPM COMPARISON
# -----------------------------

echo ">>> Comparing Salmon TPMs — EDTA"

python3 "$SCRIPT_COMPARE" \
  --blast "$EDTA_BLAST" \
  --tpm_instances "$EDTA_TPM_INST" \
  --tpm_representatives "$EDTA_TPM_REPR" \
  --min_cov 0.8 \
  --out "$EDTA_OUT_COMPARE"

echo ">>> Pipeline finished successfully"
