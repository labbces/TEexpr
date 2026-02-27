#!/usr/bin/env python3
"""
Extract genomic copies (instances) of Transposable Elements (TEs).

This script extracts genomic TE instances from a RepeatMasker .out or BED file and writes
a FASTA suitable for quantification (e.g. Salmon), along with metadata TSV and
"transcript-to-gene" map. Optionally filters by alignment coverage
(requires subject/consensus FASTA), and exports a filtered .out file.

Requirements:
    Python 3.8+
    biopython >= 1.79
    pyfaidx >= 0.7
    pandas

Install dependencies:
    pip install biopython pyranges pandas

Usage examples:
    # Basic usage with RepeatMasker output
    python extract_te_copies.py -r genome.out -g genome.fa -o te_copies.fa

    # With coverage filter (keep only instances >= 80% of consensus length)
    python extract_te_copies.py -r genome.out -g genome.fa -o te_copies.fa \\
        -s te_subject_consensi.fa -c 0.80
"""

import argparse
import sys
import re
import logging
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_repeatmasker_out(rm_path: str) -> pd.DataFrame:

    """
    Parse a RepeatMasker .out file into a DataFrame.

    Extracts query coordinates (cols 4–8), repeat name (col 9), and consensus
    alignment info (cols 11–13). Subject start/end are normalised to always be
    ascending, compensating for RepeatMasker's strand-dependent column swap.
    """

    rows = []
    with open(rm_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue

            # Skip RepeatMasker header lines
            stripped = line.lstrip()
            if stripped.startswith("SW") or stripped.startswith("score"):
                continue

            fields = line.split()

            # Hit lines start with an integer (SW score)
            if not re.match(r"^\d+$", fields[0]):
                continue

            try:
                contig  = fields[4]
                start   = int(fields[5]) - 1   # Convert to 0-based
                end     = int(fields[6])         # End is exclusive in Python slicing
                strand  = fields[8]              # '+' or 'C'
                name    = fields[9]

                # query_left is stored as "(N)" – strip the parentheses
                query_left = int(fields[7].strip("()"))

                # Subject (consensus) alignment coordinates
                # RepeatMasker swaps columns for complement hits
                s_a = int(fields[11].strip("()"))
                s_b = int(fields[12].strip("()"))
                subject_left = int(fields[13].strip("()"))

                # Normalise so aln_start <= aln_end (Using min/max, regardless of strand)
                subject_aln_start = min(s_a, s_b)
                subject_aln_end   = max(s_a, s_b)

                rows.append({
                    "Chromosome"       : contig,
                    "Start"            : start,
                    "End"              : end,
                    "Strand"           : strand,
                    "TE_name"          : name,
                    "query_left"       : query_left,
                    "subject_aln_start": subject_aln_start,
                    "subject_aln_end"  : subject_aln_end,
                    "subject_left"     : subject_left,
                })
            except Exception:
                continue   # Malformed line – skip silently

    df = pd.DataFrame(rows)

    # Normalise strand: RepeatMasker uses 'C' for complement, convert to '-'
    if not df.empty:
        df["Strand"] = df["Strand"].replace("C", "-")
    return df


# ---------------------------------------------------------------------------
# Coverage filter
# ---------------------------------------------------------------------------

def compute_coverage_filter(df, subject_fasta, min_coverage=0.80):
    """
    Filter TE instances by alignment coverage relative to consensus sequence length.

    Coverage is defined as:
        coverage = aligned_bases_in_consensus / consensus_length

    Where:
        aligned_bases_in_consensus = subject_aln_end - subject_aln_start + 1
        consensus_length           = subject_aln_end + subject_left
                                     (aligned region end + remaining bases)

    So the full consensus length can be reconstructed as:
        consensus_length = subject_aln_end + subject_left

    We use the consensus FASTA (--subject_fasta) as an independent reference for
    consensus lengths. When available, these lengths override the reconstructed
    estimate, which can be inaccurate if the annotation was trimmed.
    """

    # Build a dict of consensus lengths from the FASTA
    log.info("Loading subject/consensus FASTA for coverage estimation...")
    consensus_lengths: dict[str, int] = {}
    for rec in SeqIO.parse(subject_fasta, "fasta"):
        consensus_lengths[rec.id] = len(rec.seq)
    log.info(f"  {len(consensus_lengths)} consensus sequences loaded.")

    def _coverage(row):
        te = row["TE_name"]
        aln_start = row["subject_aln_start"]  # 1-based
        aln_end   = row["subject_aln_end"]    # 1-based, inclusive
        subj_left = row["subject_left"]

        aln_len = aln_end - aln_start + 1

        # Consensus total length: prefer FASTA lookup, fall back to reconstruction
        if te in consensus_lengths:
            cons_len = consensus_lengths[te]
        else:
            # Reconstruct: end of alignment + remaining bases
            cons_len = aln_end + subj_left
            if cons_len <= 0:
                return 0.0

        return aln_len / cons_len if cons_len > 0 else 0.0

    required_cols = {"subject_aln_start", "subject_aln_end", "subject_left"}
    missing = required_cols - set(df.columns)
    if missing:
        log.warning(
            f"Coverage filter requires RepeatMasker input. "
            f"Missing columns: {missing}. Skipping coverage filter."
        )
        return df, pd.DataFrame(columns=df.columns)

    df = df.copy()
    # TODO: replace df.apply(_coverage, axis=1) with a vectorised pandas operation to reduce runtime on large datasets.
    # Maybe use np.where() for the consensus length lookup and direct column arithmetic instead of a row-wise Python function.

    df["coverage"] = df.apply(_coverage, axis=1)

    passing = df[df["coverage"] >= min_coverage].copy()
    failing = df[df["coverage"]  < min_coverage].copy()
    return passing, failing


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
            description=(
            "Extract genomic copies of Transposable Elements from a RepeatMasker .out "
            "and write a FASTA suitable for Salmon index and quantification."),
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--repeatmasker", "-r", required=True,
        help="RepeatMasker .out annotation file.")

    parser.add_argument("--genome", "-g", required=True,
        help="Reference genome FASTA file.")

    parser.add_argument("--out_fasta", "-o", required=True,
        help=("Output FASTA file. Two additional files are created automatically: "
              "<out_fasta>.tsv (metadata) and <out_fasta>_transcript_to_gene_map.tsv."))

    parser.add_argument("--min_length", "-l", type=int, default=20,
        help="Minimum instance length in bp to retain (default: 20).")

    parser.add_argument("--subject_fasta", "-s", metavar="FASTA",
        help=("Consensus / library FASTA file (RepeatMasker subject sequences). "
              "When provided, enables the coverage filter."))

    parser.add_argument("--min_coverage", "-c", type=float, default=0.80,
        help=("Minimum alignment coverage relative to consensus length to classify an "
              "instance as 'complete' (default: 0.80 = 80%%). "
              "Only used when --subject_fasta is provided."))

    args = parser.parse_args()

    # Parse annotation file
    log.info("Parsing RepeatMasker .out...")
    df = parse_repeatmasker_out(args.repeatmasker)

    if df.empty:
        log.error("No hits found in the input file. Exiting.")
        sys.exit(1)

    log.info(f"  {len(df)} hits loaded.")

    # Minimum length filter
    df = df[(df["End"] - df["Start"]) >= args.min_length]
    log.info(f"  {len(df)} hits remaining after minimum length filter ({args.min_length} bp).")

    total_before_coverage = len(df)

    # Coverage filter (requires --subject_fasta)
    if args.subject_fasta:
        log.info(f"Applying coverage filter (min coverage: {args.min_coverage * 100:.0f}%)...")
        df_passing, df_failing = compute_coverage_filter(df, args.subject_fasta, args.min_coverage)
        log.info(
            f"  {len(df_passing)} instances passed (>= {args.min_coverage*100:.0f}% coverage); "
            f"{len(df_failing)} instances failed."
        )

        # Per-family statistics
        fam_before = df.groupby("TE_name").size().rename("before_coverage_filter")
        fam_after  = df_passing.groupby("TE_name").size().rename("after_coverage_filter")
        fam_stats  = pd.concat([fam_before, fam_after], axis=1).fillna(0).astype(int)
        fam_stats["removed"] = fam_stats["before_coverage_filter"] - fam_stats["after_coverage_filter"]
        fam_stats["pct_retained"] = (fam_stats["after_coverage_filter"] / fam_stats["before_coverage_filter"] * 100).round(1)

        stats_lines = [
            f"# Coverage filter summary (threshold: {args.min_coverage * 100:.0f}%)",
            f"# Total instances before filter : {total_before_coverage}",
            f"# Instances passing filter      : {len(df_passing)}",
            f"# Instances removed by filter   : {len(df_failing)}",
            "#",
            "\t".join(["TE_family", "before_filter", "after_filter", "removed", "pct_retained"]),
        ]
        for te_name, row_s in fam_stats.iterrows():
            stats_lines.append(
                f"{te_name}\t{row_s.before_coverage_filter}\t"
                f"{row_s.after_coverage_filter}\t{row_s.removed}\t{row_s.pct_retained}"
            )

        # Write stats file
        stats_path = args.out_fasta + ".coverage_stats.tsv"
        with open(stats_path, "w") as f:
            f.write("\n".join(stats_lines) + "\n")
        log.info(f"Coverage statistics written to: {stats_path}")

        # Write failing instances TSV (for size distribution plot and further analysis)
        failing_path = args.out_fasta + ".failing_instances.tsv"
        failing_out  = df_failing[["Chromosome", "Start", "End", "Strand", "TE_name", "coverage"]].copy()
        failing_out["length"] = failing_out["End"] - failing_out["Start"]
        failing_out.to_csv(failing_path, sep="\t", index=False)
        log.info(f"Failing instances written to: {failing_path}")

        # Write filtered RepeatMasker .out
        filtered_out_path = args.out_fasta + ".filtered.out"
        _write_filtered_rm_out(args.repeatmasker, df_passing, filtered_out_path)
        log.info(f"Filtered RepeatMasker .out written to: {filtered_out_path}")

        df = df_passing

    # Load reference genome
    log.info("Loading reference genome...")
    genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    log.info(f"  {len(genome)} sequences loaded.")

    # Extract sequences and build output
    records   = []
    tsv_lines = []
    map_lines = []
    te_counts: dict[str, int] = defaultdict(int)

    skipped = 0
    for row in df.itertuples(index=False):
        te_counts[row.TE_name] += 1
        copy_num = te_counts[row.TE_name]

        copy_id = (
            f"{row.TE_name}_copy{copy_num:04d}"
            f"|{row.Chromosome}:{row.Start + 1}-{row.End}"
            f"|{row.Strand}"
        )

        if row.Chromosome not in genome:
            log.warning(f"Contig {row.Chromosome!r} not found in genome. Skipping {copy_id}.")
            skipped += 1
            continue

        seq_obj = genome[row.Chromosome].seq[row.Start:row.End]
        seqlen  = len(seq_obj)

        if seqlen < args.min_length:
            continue

        # Reverse-complement for minus-strand hits
        if row.Strand == "-":
            seq_obj = seq_obj.reverse_complement()

        seq = str(seq_obj)

        records.append(SeqRecord(seq_obj, id=copy_id,
                description=(
                    f"family={row.TE_name} "
                    f"location={row.Chromosome}:{row.Start + 1}-{row.End}"
                    f"{row.Strand})"))
        )

        tsv_lines.append(f"{copy_id}\t{row.TE_name}\t{row.Chromosome}\t"
                         f"{row.Start + 1}\t{row.End}\t{row.Strand}\t{seqlen}")

        map_lines.append(f"{copy_id}\t{row.TE_name}")

    if skipped:
        log.warning(f"{skipped} instances skipped (contig not found in genome).")

    # Write output files
    log.info(f"Writing {len(records)} sequences to {args.out_fasta}...")
    SeqIO.write(records, args.out_fasta, "fasta")

    tsv_path = args.out_fasta + ".tsv"
    with open(tsv_path, "w") as f:
        header = "\t".join(["copy_id", "TE_name", "contig", "start", "end", "strand", "length"])
        f.write(header + "\n")
        for line in tsv_lines:
            f.write(line + "\n")
    log.info(f"Metadata TSV written to: {tsv_path}")

    map_path = args.out_fasta + "_transcript_to_gene_map.tsv"
    with open(map_path, "w") as f:
        for line in map_lines:
            f.write(line + "\n")
    log.info(f"Transcript-to-gene map written to: {map_path}")

    log.info("Done.")


# Write filtered RepeatMasker .out
def _write_filtered_rm_out(original_rm_path, passing_df, output_path):
    """
    Write a filtered RepeatMasker .out containing only instances in passing_df.
    Header lines are always preserved. Lookup is done via (contig, start, end)
    tuples, with start converted back to 1-based to match the original file.
    """
    # Build lookup set from the passing DataFrame
    # Start in df is 0-based, RepeatMasker is 1-based → Start + 1
    keep: set[tuple] = set(
        zip(
            passing_df["Chromosome"],
            passing_df["Start"] + 1,
            passing_df["End"],
        )
    )

    with open(original_rm_path) as fin, open(output_path, "w") as fout:
        for raw in fin:
            line = raw.rstrip("\n")
            if not line:
                fout.write(raw)
                continue

            stripped = line.lstrip()
            # Always keep header lines
            if stripped.startswith("SW") or stripped.startswith("score") or not stripped:
                fout.write(raw)
                continue

            fields = line.split()
            if not re.match(r"^\d+$", fields[0]):
                fout.write(raw)
                continue

            try:
                contig = fields[4]
                start  = int(fields[5])   # 1-based (original, no conversion)
                end    = int(fields[6])
                if (contig, start, end) in keep:
                    fout.write(raw)
            except (IndexError, ValueError):
                fout.write(raw)


if __name__ == "__main__":
    main()
