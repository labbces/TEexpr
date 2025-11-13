"""
RNA-Seq Processing Pipeline

This script downloads RNA-Seq data using wget, processes it with quality control,
cleans the data, and performs quantification with Salmon, then delete raw data.

The pipeline works as follows:
1. Downloads samples sequentially (one by one)
2. As soon as a sample is downloaded, it enters the parallel processing queue
3. Processing steps (QC, cleaning, quantification) run in parallel with N workers
4. Delete raw RNA-Seq data

Usage:
    python pipeline_init.py --sra_list samples.txt --transcripts transcriptome.fa --outdir output --max_workers 4
"""

import os
import subprocess
import argparse
import glob
import multiprocessing as mp
import threading
import queue
import time
import sys
import shutil
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--sra_list", required=True, help="File with the IDs (SRR)")
parser.add_argument("--outdir", default="pipeline_output", help="Output directory")
parser.add_argument("--max_parallel", default="4", help="Number of parallel sample processing for each tools")
parser.add_argument("--transcripts", required=True, help="Path to transcriptome FASTA file")
parser.add_argument("--salmon_threads", default="10", help="Number of threads for Salmon")
parser.add_argument("--max_workers", default="4", help="Maximum number of parallel workers")
parser.add_argument("--base_url", default="https://your-data-source.com/fastq", 
                    help="Base URL for downloading FASTQ files")
args = parser.parse_args()

# Creating directories
outdir = os.path.abspath(os.path.expanduser(args.outdir))
raw_dir = os.path.join(outdir, "raw_data")
fastqc_dir = os.path.join(outdir, "fastqc_results")
clean_dir = os.path.join(outdir, "clean_data")
salmon_index_dir = os.path.join(outdir, "salmon_index")
salmon_results_dir = os.path.join(outdir, "salmon_results")

temp_dir = "/tmp"

os.makedirs(outdir, exist_ok=True)
os.makedirs(raw_dir, exist_ok=True)
os.makedirs(fastqc_dir, exist_ok=True)
os.makedirs(clean_dir, exist_ok=True)
os.makedirs(salmon_index_dir, exist_ok=True)
os.makedirs(salmon_results_dir, exist_ok=True)

# Read IDs and fill list with SRRs
sra_ids = []
with open(args.sra_list) as f:
    for line in f:
        sra_ids.append(line.strip())

# Global variables for tracking progress
download_queue = queue.Queue()
processed_samples = {}
salmon_index_built = threading.Event()
processing_stats = {
    'downloaded': 0,
    'fastqc_completed': 0,
    'cleaned': 0,
    'quantified': 0,
    'failed_downloads': 0,
    'failed_processing': 0
}
stats_lock = threading.Lock()

def update_stats(category, increment=1):
    """Thread-safe update of processing statistics"""
    with stats_lock:
        processing_stats[category] += increment

def check_download_completed(srr):
    """Check if download is already completed for a sample"""
    r1_file = os.path.join(raw_dir, f"{srr}_1.fastq")
    r2_file = os.path.join(raw_dir, f"{srr}_2.fastq")
    return os.path.exists(r1_file) and os.path.exists(r2_file)

def check_fastqc_completed(srr):
    """Check if FastQC is already completed for a sample"""
    r1_html = os.path.join(fastqc_dir, f"{srr}_1_fastqc.html")
    r2_html = os.path.join(fastqc_dir, f"{srr}_2_fastqc.html")
    return os.path.exists(r1_html) and os.path.exists(r2_html)

def check_cleaning_completed(srr):
    """Check if BBduk cleaning is already completed for a sample"""
    clean_r1 = os.path.join(clean_dir, f"{srr}_clean_1.fastq")
    clean_r2 = os.path.join(clean_dir, f"{srr}_clean_2.fastq")
    return os.path.exists(clean_r1) and os.path.exists(clean_r2)

def check_quantification_completed(srr):
    """Check if Salmon quantification is already completed for a sample"""
    quant_file = os.path.join(salmon_results_dir, srr, "quant.sf")
    return os.path.exists(quant_file)

def download_sample(srr):
    """Download a single sample using wget"""
    # Check if download already completed
    if check_download_completed(srr):
        print(f"[DOWNLOAD] ✓ {srr} already downloaded, skipping...")
        update_stats('downloaded')
        return True
    
    print(f"[DOWNLOAD] Starting download of {srr}...")
    
    # Create directory for this sample
    sample_dir = os.path.join(raw_dir, srr)
    os.makedirs(sample_dir, exist_ok=True)
    
    try:
        # Download R1 and R2 files using wget
        r1_url = f"{args.base_url}/{srr}_1.fastq.gz"
        r2_url = f"{args.base_url}/{srr}_2.fastq.gz"
        
        r1_output = os.path.join(sample_dir, f"{srr}_1.fastq.gz")
        r2_output = os.path.join(sample_dir, f"{srr}_2.fastq.gz")
        
        # Download R1
        wget_cmd_r1 = [
            "wget", "-O", r1_output, "--progress=bar:force",
            "--timeout=30", "--tries=3", r1_url
        ]
        
        result_r1 = subprocess.run(wget_cmd_r1, capture_output=True, text=True)
        if result_r1.returncode != 0:
            print(f"[DOWNLOAD] Error downloading R1 for {srr}: {result_r1.stderr}")
            return False
        
        # Download R2
        wget_cmd_r2 = [
            "wget", "-O", r2_output, "--progress=bar:force",
            "--timeout=30", "--tries=3", r2_url
        ]
        
        result_r2 = subprocess.run(wget_cmd_r2, capture_output=True, text=True)
        if result_r2.returncode != 0:
            print(f"[DOWNLOAD] Error downloading R2 for {srr}: {result_r2.stderr}")
            if os.path.exists(r1_output):
                os.remove(r1_output)
            return False
        
        # Decompress files
        decomp_r1 = os.path.join(sample_dir, f"{srr}_1.fastq")
        decomp_r2 = os.path.join(sample_dir, f"{srr}_2.fastq")
        
        # Decompress R1
        with open(decomp_r1, 'w') as f:
            gunzip_result_r1 = subprocess.run(["gunzip", "-c", r1_output], 
                                            stdout=f, stderr=subprocess.PIPE, text=True)
        if gunzip_result_r1.returncode != 0:
            print(f"[DOWNLOAD] Error decompressing R1 for {srr}")
            return False
        
        # Decompress R2
        with open(decomp_r2, 'w') as f:
            gunzip_result_r2 = subprocess.run(["gunzip", "-c", r2_output], 
                                            stdout=f, stderr=subprocess.PIPE, text=True)
        if gunzip_result_r2.returncode != 0:
            print(f"[DOWNLOAD] Error decompressing R2 for {srr}")
            return False
        
        # Remove compressed files and move to final location
        os.remove(r1_output)
        os.remove(r2_output)
        
        final_r1 = os.path.join(raw_dir, f"{srr}_1.fastq")
        final_r2 = os.path.join(raw_dir, f"{srr}_2.fastq")
        
        os.rename(decomp_r1, final_r1)
        os.rename(decomp_r2, final_r2)
        os.rmdir(sample_dir)
        
        print(f"[DOWNLOAD] ✓ Successfully downloaded {srr}")
        update_stats('downloaded')
        return True
        
    except Exception as e:
        print(f"[DOWNLOAD] ✗ Failed to download {srr}: {str(e)}")
        update_stats('failed_downloads')
        return False

def process_sample(srr):
    """Process a single sample: FastQC -> BBduk -> Salmon quantification"""
    try:
        # Check if all processing is already completed
        if (check_fastqc_completed(srr) and 
            check_cleaning_completed(srr) and 
            check_quantification_completed(srr)):
            print(f"[PROCESSING] ✓ {srr} already fully processed, skipping all steps...")
            update_stats('fastqc_completed')
            update_stats('cleaned')
            update_stats('quantified')
            return True
        
        # Step 1: FastQC
        if check_fastqc_completed(srr):
            print(f"[FASTQC] ✓ {srr} FastQC already completed, skipping...")
            update_stats('fastqc_completed')
        else:
            print(f"[FASTQC] Starting FastQC for {srr}...")
            
            r1_file = os.path.join(raw_dir, f"{srr}_1.fastq")
            r2_file = os.path.join(raw_dir, f"{srr}_2.fastq")
            
            if not (os.path.exists(r1_file) and os.path.exists(r2_file)):
                print(f"[FASTQC] ✗ Raw files not found for {srr}")
                update_stats('failed_processing')
                return False
            
            # Run FastQC
            for fastq_file in [r1_file, r2_file]:
                fastqc_cmd = ["fastqc", fastq_file, "-o", fastqc_dir, "--quiet"]
                fastqc_result = subprocess.run(fastqc_cmd, capture_output=True, text=True)
                if fastqc_result.returncode != 0:
                    print(f"[FASTQC] Warning: FastQC failed for {os.path.basename(fastq_file)}")
            
            print(f"[FASTQC] ✓ FastQC completed for {srr}")
            update_stats('fastqc_completed')
        
        # Step 2: BBduk cleaning
        if check_cleaning_completed(srr):
            print(f"[BBDUK] ✓ {srr} cleaning already completed, skipping...")
            update_stats('cleaned')
        else:
            print(f"[BBDUK] Starting cleaning for {srr}...")
            
            r1_file = os.path.join(raw_dir, f"{srr}_1.fastq")
            r2_file = os.path.join(raw_dir, f"{srr}_2.fastq")
            
            if not (os.path.exists(r1_file) and os.path.exists(r2_file)):
                print(f"[BBDUK] ✗ Raw files not found for {srr}")
                update_stats('failed_processing')
                return False
            
            clean_r1 = os.path.join(clean_dir, f"{srr}_clean_1.fastq")
            clean_r2 = os.path.join(clean_dir, f"{srr}_clean_2.fastq")
            
            bbduk_cmd = [
                "bbduk.sh", f"in1={r1_file}", f"in2={r2_file}",
                f"out1={clean_r1}", f"out2={clean_r2}",
                "ref=adapters", "ktrim=r", "k=23", "mink=11", "hdist=1",
                "tpe", "tbo", "qtrim=rl", "trimq=10", "minlen=50",
                "threads=" + args.threads
            ]
            
            bbduk_result = subprocess.run(bbduk_cmd, capture_output=True, text=True)
            if bbduk_result.returncode != 0:
                print(f"[BBDUK] ✗ Error cleaning {srr}: {bbduk_result.stderr}")
                update_stats('failed_processing')
                return False
            
            print(f"[BBDUK] ✓ Successfully cleaned {srr}")
            update_stats('cleaned')
        
        # Step 3: Salmon quantification
        if check_quantification_completed(srr):
            print(f"[SALMON] ✓ {srr} quantification already completed, skipping...")
            update_stats('quantified')
        else:
            # Wait for index to be built
            print(f"[SALMON] Waiting for index to be built...")
            salmon_index_built.wait()  # Wait for index to be ready
            
            print(f"[SALMON] Starting quantification for {srr}...")
            
            clean_r1 = os.path.join(clean_dir, f"{srr}_clean_1.fastq")
            clean_r2 = os.path.join(clean_dir, f"{srr}_clean_2.fastq")
            
            if not (os.path.exists(clean_r1) and os.path.exists(clean_r2)):
                print(f"[SALMON] ✗ Clean files not found for {srr}")
                update_stats('failed_processing')
                return False
            
            sample_output = os.path.join(salmon_results_dir, srr)
            
            salmon_quant_cmd = [
                "salmon", "quant", "-i", salmon_index_dir, "-l", "A",
                "-1", clean_r1, "-2", clean_r2,
                "-p", args.salmon_threads, "-o", sample_output
            ]
            
            quant_result = subprocess.run(salmon_quant_cmd, capture_output=True, text=True)
            if quant_result.returncode != 0:
                print(f"[SALMON] ✗ Error quantifying {srr}: {quant_result.stderr}")
                update_stats('failed_processing')
                return False
            
            print(f"[SALMON] ✓ Successfully quantified {srr}")
            update_stats('quantified')
        
        # Clean up raw files for this sample to save space (only if not already cleaned)
        r1_file = os.path.join(raw_dir, f"{srr}_1.fastq")
        r2_file = os.path.join(raw_dir, f"{srr}_2.fastq")
        
        if os.path.exists(r1_file) or os.path.exists(r2_file):
            try:
                if os.path.exists(r1_file):
                    os.remove(r1_file)
                if os.path.exists(r2_file):
                    os.remove(r2_file)
                print(f"[CLEANUP] ✓ Removed raw files for {srr}")
            except Exception as e:
                print(f"[CLEANUP] Warning: Could not remove raw files for {srr}: {e}")
        
        return True
        
    except Exception as e:
        print(f"[PROCESSING] ✗ Failed processing {srr}: {str(e)}")
        update_stats('failed_processing')
        return False

def build_salmon_index():
    """Build Salmon index in a separate thread"""
    # Check if index already exists
    index_info_file = os.path.join(salmon_index_dir, "info.json")
    if os.path.exists(index_info_file):
        print("[INDEX] ✓ Salmon index already exists, skipping build...")
        salmon_index_built.set()  # Signal that index is ready
        return
    
    print("[INDEX] Building Salmon index...")
    
    if not os.path.exists(args.transcripts):
        print(f"[INDEX] ✗ Transcriptome file not found: {args.transcripts}")
        sys.exit(1)
    
    salmon_index_cmd = [
        "salmon", "index", "-t", args.transcripts, "-i", salmon_index_dir
    ]
    
    index_result = subprocess.run(salmon_index_cmd, capture_output=True, text=True)
    if index_result.returncode != 0:
        print(f"[INDEX] ✗ Error building Salmon index: {index_result.stderr}")
        sys.exit(1)
    
    print("[INDEX] ✓ Salmon index successfully created!")
    salmon_index_built.set()  # Signal that index is ready

# Main pipeline execution
print("=" * 80)
print("PARALLEL RNA-SEQ PROCESSING PIPELINE")
print("=" * 80)
print(f"Total samples: {len(sra_ids)}")
print(f"Max parallel workers: {args.max_workers}")
print(f"Base URL: {args.base_url}")

# Check for existing results and show resume information
print("\n" + "=" * 80)
print("CHECKING EXISTING RESULTS")
print("=" * 80)

already_downloaded = 0
already_fastqc = 0
already_cleaned = 0
already_quantified = 0

for srr in sra_ids:
    if check_download_completed(srr):
        already_downloaded += 1
    if check_fastqc_completed(srr):
        already_fastqc += 1
    if check_cleaning_completed(srr):
        already_cleaned += 1
    if check_quantification_completed(srr):
        already_quantified += 1

print(f"Found existing results:")
print(f"  - Downloaded: {already_downloaded}/{len(sra_ids)} samples")
print(f"  - FastQC completed: {already_fastqc}/{len(sra_ids)} samples")
print(f"  - Cleaned: {already_cleaned}/{len(sra_ids)} samples")
print(f"  - Quantified: {already_quantified}/{len(sra_ids)} samples")

if already_quantified == len(sra_ids):
    print(f"\n✓ All samples already fully processed! Skipping to merge step...")
    # Update stats with existing results
    processing_stats['downloaded'] = already_downloaded
    processing_stats['fastqc_completed'] = already_fastqc
    processing_stats['cleaned'] = already_cleaned
    processing_stats['quantified'] = already_quantified
    # Set index as built (assume it exists if quantification is done)
    salmon_index_built.set()
else:
    print(f"\n→ Will process {len(sra_ids) - already_quantified} remaining samples")

print("=" * 80)

# Start building Salmon index in background
index_thread = threading.Thread(target=build_salmon_index)
index_thread.daemon = True
index_thread.start()

# Start the parallel processing pipeline only if needed
if already_quantified < len(sra_ids):
    with ThreadPoolExecutor(max_workers=int(args.max_workers)) as executor:
        futures = []
        
        for i, srr in enumerate(sra_ids, 1):
            # Skip if already fully processed
            if check_quantification_completed(srr):
                print(f"[PIPELINE] ✓ {srr} already fully processed, skipping...")
                continue
                
            # Download sequentially (only if not already downloaded)
            print(f"\n[PIPELINE] Processing sample {i}/{len(sra_ids)}: {srr}")
            
            # Download the sample
            if download_sample(srr):
                # Submit for parallel processing
                future = executor.submit(process_sample, srr)
                futures.append((srr, future))
                print(f"[PIPELINE] ✓ {srr} submitted for parallel processing")
            else:
                print(f"[PIPELINE] ✗ {srr} download failed, skipping processing")
        
        # Wait for all processing to complete
        if futures:
            print(f"\n[PIPELINE] Waiting for all {len(futures)} submitted samples to complete processing...")
            
            for srr, future in futures:
                try:
                    result = future.result()  # This will block until the future completes
                    if result:
                        print(f"[PIPELINE] ✓ {srr} processing completed successfully")
                    else:
                        print(f"[PIPELINE] ✗ {srr} processing failed")
                except Exception as e:
                    print(f"[PIPELINE] ✗ {srr} processing exception: {str(e)}")
        else:
            print("\n[PIPELINE] ✓ No samples needed processing, all already completed!")
else:
    print("\n[PIPELINE] ✓ All samples already processed, skipping pipeline execution!")

# Wait for index building to complete
index_thread.join()

print("\n" + "=" * 80)
print("MERGING QUANTIFICATION RESULTS")
print("=" * 80)

# Merge quantification results
if processing_stats['quantified'] > 0:
    print("Merging quantification results...")
    
    # Find all directories with quant.sf files
    quant_dirs = []
    for item in os.listdir(salmon_results_dir):
        quant_file = os.path.join(salmon_results_dir, item, "quant.sf")
        if os.path.isfile(quant_file):
            quant_dirs.append(os.path.join(salmon_results_dir, item))
    
    if quant_dirs:
        # Merge TPM values
        tpm_merge_cmd = [
            "salmon", "quantmerge",
            "--column", "TPM",
            "--output", os.path.join(outdir, "merged_TPM.tsv"),
            "--quants"
        ] + quant_dirs
        
        tpm_result = subprocess.run(tpm_merge_cmd, capture_output=True, text=True)
        if tpm_result.returncode != 0:
            print(f"[MERGE] ✗ Error merging TPM: {tpm_result.stderr}")
        else:
            print("[MERGE] ✓ TPM merge completed successfully")
        
        # Merge count values
        counts_merge_cmd = [
            "salmon", "quantmerge",
            "--column", "NumReads",
            "--output", os.path.join(outdir, "merged_counts.tsv"),
            "--quants"
        ] + quant_dirs
        
        counts_result = subprocess.run(counts_merge_cmd, capture_output=True, text=True)
        if counts_result.returncode != 0:
            print(f"[MERGE] ✗ Error merging counts: {counts_result.stderr}")
        else:
            print("[MERGE] ✓ Counts merge completed successfully")

print("\n" + "=" * 80)
print("FINAL CLEANUP")
print("=" * 80)

# Clean up remaining raw and clean directories
print("[CLEANUP] Removing remaining temporary files...")
try:
    if os.path.exists(clean_dir):
        clean_size = sum(os.path.getsize(os.path.join(dirpath, filename))
                        for dirpath, dirnames, filenames in os.walk(clean_dir)
                        for filename in filenames) / (1024*1024)  # Size in MB
        
        shutil.rmtree(clean_dir)
        print(f"[CLEANUP] ✓ Removed clean data directory ({clean_size:.1f} MB freed)")
        
except Exception as e:
    print(f"[CLEANUP] Warning: Error during cleanup: {e}")

print("\n" + "=" * 80)
print("PIPELINE COMPLETED!")
print("=" * 80)

# Final statistics
print(f"Processing Statistics:")
print(f"  Downloaded: {processing_stats['downloaded']}/{len(sra_ids)}")
print(f"  FastQC completed: {processing_stats['fastqc_completed']}/{len(sra_ids)}")
print(f"  Cleaned: {processing_stats['cleaned']}/{len(sra_ids)}")
print(f"  Quantified: {processing_stats['quantified']}/{len(sra_ids)}")
print(f"  Failed downloads: {processing_stats['failed_downloads']}")
print(f"  Failed processing: {processing_stats['failed_processing']}")

print(f"\nResults directory: {outdir}")
if processing_stats['downloaded'] > 0:
    print(f"  - FastQC results: {fastqc_dir}")
    print(f"  - Salmon quantification results: {salmon_results_dir}")
    if os.path.exists(os.path.join(outdir, "merged_TPM.tsv")):
        print(f"  - Merged TPM results: {os.path.join(outdir, 'merged_TPM.tsv')}")
    if os.path.exists(os.path.join(outdir, "merged_counts.tsv")):
        print(f"  - Merged counts results: {os.path.join(outdir, 'merged_counts.tsv')}")

# Exit with appropriate code
total_failures = processing_stats['failed_downloads'] + processing_stats['failed_processing']
sys.exit(0 if total_failures == 0 else 1)
