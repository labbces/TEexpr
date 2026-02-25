#!/usr/bin/env python3
"""
Estimate sequence similarity between RepeatMasker genomic instances (queries)
and their TE representative (consensus) using base-pair level calculation.

Identity is calculated as:
    matches = alignment_length - mismatches - indels
    identity = matches / alignment_length * 100

Bins are automatically defined from percentiles of the identity distribution.
Sequences with negative identity are saved separately for manual inspection.
SW score is included as a column for reference or filtering.
"""

import argparse
import sys
import re
import numpy as np

# Parse RepeatMasker .out file and calculate identity per instance
def parse_repeatmasker_out(rm_out, negative_file):
    identities = []
    results = []

    with open(rm_out) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            # Skip headers
            if line.lstrip().startswith("SW") or line.lstrip().startswith("score"):
                continue
            fields = line.split()
            # Valid RepeatMasker hits start with a numeric SW score
            if not re.match(r"^\d+$", fields[0]):
                continue

            try:
                sw_score = float(fields[0])
                perc_div = float(fields[1])
                perc_del = float(fields[2])
                perc_ins = float(fields[3])

                chrom = fields[4]
                start = int(fields[5])
                end = int(fields[6])
                strand = fields[8]
                te_name = fields[9]

                copy_id = f"{te_name}|{chrom}:{start}-{end}|{strand}"

                # Alignment length
                alignment_length = end - start + 1

                # Convert percentages to counts of bases
                mismatches_bp = perc_div / 100 * alignment_length
                indels_bp = (perc_ins + perc_del) / 100 * alignment_length
                matches_bp = alignment_length - mismatches_bp - indels_bp

                # Raw identity before sanity checks
                identity_raw = (matches_bp / alignment_length) * 100

                # Handle negative identity
                if identity_raw < 0:
                    print(f"Warning (line {lineno}): negative identity ({identity_raw:.2f}) for {copy_id}", file=sys.stderr)
                    negative_file.write("\t".join([
                        copy_id, te_name, f"{identity_raw:.2f}",
                        str(alignment_length), str(perc_div), str(perc_ins), str(perc_del),
                        str(sw_score)
                    ]) + "\n")
                    identity = 0.0
                elif identity_raw > 100:
                    identity = 100.0
                else:
                    identity = identity_raw

                results.append((copy_id, te_name, round(identity, 2), sw_score))
                identities.append(identity)

            except Exception as e:
                print(f"Warning (line {lineno}): could not parse line:\n{line}\nReason: {e}", file=sys.stderr)
                continue

    return results, np.array(identities)

# Assign identity to automatic similarity bins using percentiles
def assign_bins(values, n_bins):
    # Compute percentile edges
    percentiles = np.percentile(values, np.linspace(0, 100, n_bins + 1))
    bins = []
    for i in range(n_bins):
        bins.append((percentiles[i], percentiles[i + 1]))
    return bins

def main():
    parser = argparse.ArgumentParser(description="Estimate similarity between RepeatMasker instances and TE representatives")
    parser.add_argument("-i", "--input", required=True, help="RepeatMasker .out file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument("-b", "--bins", type=int, default=4, help="Number of bins based on percentiles (default 4)")
    args = parser.parse_args()

    # Open file to store sequences with negative identity
    negative_file_path = args.output + ".negative.txt"
    negative_file = open(negative_file_path, "w")
    negative_file.write("\t".join(["copy_id","TE_name","identity_raw","alignment_length","perc_div","perc_ins","perc_del","sw_score"]) + "\n")

    # Parse RepeatMasker output
    results, identities = parse_repeatmasker_out(args.input, negative_file)
    negative_file.close()

    # Compute automatic bins
    bins = assign_bins(identities, args.bins)

    # Assign each instance to a bin
    final_results = []
    for copy_id, te_name, identity, sw_score in results:
        for low, high in bins:
            if low <= identity <= high or (high == 100 and identity <= 100):
                bin_label = f"{int(round(low))}-{int(round(high))}"
                break
        final_results.append((copy_id, te_name, identity, bin_label, sw_score))

    # Write output
    with open(args.output, "w") as out:
        out.write("\t".join(["copy_id", "TE_name", "estimated_identity", "similarity_bin", "sw_score"]) + "\n")
        for r in final_results:
            out.write("\t".join(map(str, r)) + "\n")

    print(f"Finished. {len(final_results)} entries written to {args.output}")
    print(f"Negative identities saved in {negative_file_path}")

if __name__ == "__main__":
    main()
