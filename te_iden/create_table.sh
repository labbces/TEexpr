#!/bin/bash/env python3
#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 4
#$ -l h=neotera

module load python3

# Input directories
FLTE_DIR="/Storage/data2/andreza.cunha/results/flTE/"
DEEPTE_DIR="/Storage/data2/andreza.cunha/results/DeepTE"
TESORTER_DIR="/Storage/data2/andreza.cunha/results/TE_sorter"

# Output directory
OUTPUT_DIR="/Storage/data2/andreza.cunha/results/te_overview_table"
mkdir -p "$OUTPUT_DIR"

# Combined output file
COMBINED_FILE="$OUTPUT_DIR/all_species_combined.csv"
> "$COMBINED_FILE"

# Species list (prefixes)
SPECIES=("GreenFoxtail" "Shybrid" "Sofficinarum" "Sspontaneum" "cane19" "cane87" "caneNp" "canek3" "maize" "sorghum")

# Path to Python script
PYTHON_SCRIPT="/Storage/data2/andreza.cunha/scripts/classification_table.py"

# Add header to combined file
echo "TE_ID,EarlGrey,DeepTE,TEsorter,Species" > "$COMBINED_FILE"

# Loop over each species
for SPECIES_NAME in "${SPECIES[@]}"; do
    echo "Processing $SPECIES_NAME..."

    # Input files
    FLTE_FILE="$FLTE_DIR/${SPECIES_NAME}.flTE.fa"
    DEEPTE_FILE="$DEEPTE_DIR/${SPECIES_NAME}_opt_DeepTE.fasta"
    TESORTER_FILE="$TESORTER_DIR/${SPECIES_NAME}.cls.tsv"

    # Check if files exist
    if [[ -f "$FLTE_FILE" && -f "$DEEPTE_FILE" && -f "$TESORTER_FILE" ]]; then
        # Individual output file
        OUTPUT_FILE="$OUTPUT_DIR/${SPECIES_NAME}_output.csv"

        # Execute Python script
        python3 "$PYTHON_SCRIPT" \
            --earlgrey "$FLTE_FILE" \
            --deepte "$DEEPTE_FILE" \
            --tesorter "$TESORTER_FILE" \
            --output "$OUTPUT_FILE"

        echo "Result saved in $OUTPUT_FILE"

        # Add data to combined file (ignore header from individual files)
        tail -n +2 "$OUTPUT_FILE" | dos2unix | awk -v species="$SPECIES_NAME" -F',' 'BEGIN {OFS=","} {print $1, $2, $3, "\"" $4 "\"", species}' >> "$COMBINED_FILE"
    else
        echo "Files missing for $SPECIES_NAME. Skipping..."
    fi
done

# Remove ^M characters from combined file (force Unix format)
sed -i 's/\r$//' "$COMBINED_FILE"

# Sort combined file by TE_ID column (first column)
tail -n +2 "$COMBINED_FILE" | sort -t, -k1,1 > "${COMBINED_FILE}.sorted"

# Reinsert header in sorted file
head -n 1 "$COMBINED_FILE" > "${COMBINED_FILE}.tmp"
cat "${COMBINED_FILE}.sorted" >> "${COMBINED_FILE}.tmp"
mv "${COMBINED_FILE}.tmp" "$COMBINED_FILE"
rm "${COMBINED_FILE}.sorted"

echo "Combined table generated in $COMBINED_FILE"

