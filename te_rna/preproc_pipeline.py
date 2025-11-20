#!/usr/bin/python

"""
RNA-Seq Processing Pipeline

This script downloads RNA-Seq data using wget, processes it with quality control,
cleans the data, performs quantification with Salmon and delete raw data.

The pipeline works as follows:
1. Downloads samples sequentially (one by one)
2. As soon as a sample is downloaded, it enters the parallel processing queue
3. Processing steps (QC, cleaning, quantification) run in parallel with N workers

Usage:
	python pipeline_init.py --sra_list samples.txt --transcripts transcriptome.fa --outdir output --parallel_jobs 4
"""

import os
import subprocess
import argparse
import threading
import queue
import time
import sys
import shutil
import logging

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--sra_list", required=True, help="File with the IDs (SRR)")
parser.add_argument("--transcripts", required=True, help="Path to transcriptome FASTA file")
parser.add_argument("--config", default="config.yaml", help="Path to configuration YAML file")
parser.add_argument("--outdir", default="pipeline_output", help="Output directory")
parser.add_argument("--salmon_path", default="salmon", help="Path to Salmon")
parser.add_argument("--fastqc_path", default="fastqc", help="Path to FastQC")
parser.add_argument("--bbduk_path", default="bbduk.sh", help="Path to BBduk")
parser.add_argument("--parallel_jobs", default="5", help="Maximum number of parallel jobs")
parser.add_argument("--decoy", help="Path to decoy sequences FASTA file (optional)")
args = parser.parse_args()

if not os.path.exists(args.config):
	print(f"Configuration file not found: {args.config}")
	sys.exit(1)

# Parse simple config file (base_url, user, password)
with open(args.config, 'r') as config_file:
	args.base_url = None
	args.user = None
	args.password = None
	
	for line in config_file:
		line = line.strip()
		if line and not line.startswith('#'):
			if line.startswith('base_url:'):
				args.base_url = line.split(':', 1)[1].strip().strip('"\'')
			elif line.startswith('user:'):
				args.user = line.split(':', 1)[1].strip().strip('"\'')
			elif line.startswith('password:'):
				args.password = line.split(':', 1)[1].strip().strip('"\'')

# Setup logging
def setup_logging(outdir):
	"""Configure logging for the pipeline"""
	log_dir = os.path.join(outdir, "logs")
	os.makedirs(log_dir, exist_ok=True)
	
	# Create formatters
	detailed_formatter = logging.Formatter(
		'%(asctime)s - %(name)s - %(levelname)s - [%(threadName)s] - %(message)s',
		datefmt='%Y-%m-%d %H:%M:%S'
	)
	simple_formatter = logging.Formatter(
		'%(asctime)s - %(levelname)s - %(message)s',
		datefmt='%H:%M:%S'
	)
	
	# Setup root logger
	logger = logging.getLogger()
	logger.setLevel(logging.INFO)
	
	# File handler for detailed logs
	file_handler = logging.FileHandler(os.path.join(log_dir, "pipeline.log"))
	file_handler.setLevel(logging.INFO)
	file_handler.setFormatter(detailed_formatter)
	logger.addHandler(file_handler)
	
	# Console handler for user-friendly output
	console_handler = logging.StreamHandler(sys.stdout)
	console_handler.setLevel(logging.INFO)
	console_handler.setFormatter(simple_formatter)
	logger.addHandler(console_handler)
	
	# Error log file
	error_handler = logging.FileHandler(os.path.join(log_dir, "errors.log"))
	error_handler.setLevel(logging.ERROR)
	error_handler.setFormatter(detailed_formatter)
	logger.addHandler(error_handler)
	
	return logger

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

# Initialize logging
logger = setup_logging(outdir)
logger.info("Pipeline initialized successfully")
logger.info(f"Output directory: {outdir}")

# Read IDs and fill list with SRRs
sra_ids = {}
with open(args.sra_list) as f:
	for line in f:
		parts = line.strip().split()
		if parts[0] not in sra_ids:
			sra_ids[parts[0]] = [parts[1]]
		else:
			sra_ids[parts[0]].append(parts[1])

# Global variables for tracking progress
sample_queue = queue.Queue()
fastqc_queue = queue.Queue()
bbduk_queue = queue.Queue()
salmon_queue = queue.Queue()
salmon_index_built = threading.Event()
pipeline_shutdown = threading.Event()

processing_stats = {
	'downloaded': 0,
	'fastqc_completed': 0,
	'cleaned': 0,
	'quantified': 0,
	'failed_downloads': 0,
	'failed_processing': 0
}
stats_lock = threading.Lock()

# Worker counters
active_fastqc_workers = 0
active_bbduk_workers = 0
active_salmon_workers = 0
worker_lock = threading.Lock()

def update_stats(category, increment=1):
	"""Thread-safe update of processing statistics"""
	with stats_lock:
		processing_stats[category] += increment

def check_fastqc_completed(srr):
	"""Check if FastQC is already completed for a sample"""
	r1_html = os.path.join(fastqc_dir, f"{srr}_1_fastqc.html")
	r2_html = os.path.join(fastqc_dir, f"{srr}_2_fastqc.html")
	return os.path.exists(r1_html) and os.path.getsize(r1_html) > 0 and os.path.exists(r2_html) and os.path.getsize(r2_html) > 0

def fastqc_passed_quality(srr):
	"""Check if FastQC quality is sufficient (>80% good sequences) from HTML report"""
	import re
	r1_html = os.path.join(fastqc_dir, f"{srr}_1_fastqc.html")
	r2_html = os.path.join(fastqc_dir, f"{srr}_2_fastqc.html")
	def parse_html(file):
		total = None
		poor = None
		if not os.path.exists(file):
			return None, None
		import re
		for line in open(file):
			# Exemplo: <tr><td>Total Sequences</td><td>23400214</td></tr>
			m_total = re.search(r'<tr><td>Total Sequences</td><td>(\d+)</td></tr>', line)
			if m_total:
				total = int(m_total.group(1))
			m_poor = re.search(r'<tr><td>Sequences flagged as poor quality</td><td>(\d+)</td></tr>', line)
			if m_poor:
				poor = int(m_poor.group(1))
		return total, poor
	t1, p1 = parse_html(r1_html)
	t2, p2 = parse_html(r2_html)
	if t1 is None or t2 is None or p1 is None or p2 is None:
		return False
	good1 = (t1 - p1) / t1 if t1 > 0 else 0
	good2 = (t2 - p2) / t2 if t2 > 0 else 0
	return good1 > 0.8 and good2 > 0.8

def check_cleaning_completed(srr):
	"""Check if BBduk cleaning is already completed for a sample"""
	clean_r1 = os.path.join(clean_dir, f"{srr}_clean_1.fastq.gz")
	clean_r2 = os.path.join(clean_dir, f"{srr}_clean_2.fastq.gz")
	return os.path.exists(clean_r1) and os.path.getsize(clean_r1) > 0 and os.path.exists(clean_r2) and os.path.getsize(clean_r2) > 0

def check_quantification_completed(srr):
	"""Check if Salmon quantification is already completed for a sample"""
	quant_log = os.path.join(salmon_results_dir, srr, "logs", "salmon_quant.log")
	quant_file = os.path.join(salmon_results_dir, srr, "quant.sf")

	if not os.path.exists(quant_log):
		return False

	with open(quant_log, 'r') as f:
		log_content = f.read()
		if "writing output" not in log_content:
			return False

	return os.path.exists(quant_file) and os.path.getsize(quant_file) > 0

def download_sample(srr):
	"""Download a single sample using wget"""
	logger = logging.getLogger(__name__)
	logger.info(f"Starting download of {srr}")
	
	try:
		# Download R1 and R2 files using wget
		r1_url = f"{args.base_url}{srr}_1.fastq.gz"
		r2_url = f"{args.base_url}{srr}_2.fastq.gz"
		
		r1_output = os.path.join(raw_dir, f"{srr}_1.fastq.gz")
		r2_output = os.path.join(raw_dir, f"{srr}_2.fastq.gz")
		
		# Download R1
		wget_cmd_r1 = [
			"wget", "-c", "-O", r1_output, "--user", args.user, "--password", args.password,
			"--progress=bar:force", "--timeout=600", "--tries=3", r1_url
		]
		
		result_r1 = subprocess.run(wget_cmd_r1, capture_output=True, text=True)
			
		if result_r1.returncode != 0:
			if result_r1.returncode == 416:  # 416 means file already fully downloaded
				logger.info(f"{srr} R1 already downloaded, skipping")
				pass
			else:
				logger.error(f"Error downloading R1 for {srr}: {result_r1.stderr}")
				return False
		
		# Download R2
		wget_cmd_r2 = [
			"wget", "-c", "-O", r2_output, "--user", args.user, "--password", args.password,
			"--progress=bar:force", "--timeout=600", "--tries=3", r2_url
		]
		
		result_r2 = subprocess.run(wget_cmd_r2, capture_output=True, text=True)
		if result_r2.returncode != 0:
			if result_r2.returncode == 416:  # 416 means file already fully downloaded
				logger.info(f"{srr} R2 already downloaded, skipping")
				pass
			else:
				logger.error(f"Error downloading R2 for {srr}: {result_r2.stderr}")
				if os.path.exists(r1_output):
					os.remove(r1_output)
				return False
		
		logger.info(f"Successfully downloaded {srr}")
		update_stats('downloaded')
		return True
		
	except Exception as e:
		logger.error(f"Failed to download {srr}: {str(e)}")
		update_stats('failed_downloads')
		return False

def fastqc_worker():
	"""Worker for FastQC processing"""
	global active_fastqc_workers
	logger = logging.getLogger(__name__)
	
	while not pipeline_shutdown.is_set():
		try:
			species, srr = fastqc_queue.get(timeout=1)
			if srr is None:  # Shutdown signal
				break
				
			with worker_lock:
				active_fastqc_workers += 1
				
			if check_fastqc_completed(srr):
				logger.info(f"FastQC already completed for {srr}, skipping")
				update_stats('fastqc_completed')
				bbduk_queue.put([species, srr])  # Move to next stage
			else:
				logger.info(f"Starting FastQC for {srr}")
				
				r1_file = os.path.join(raw_dir, f"{srr}_1.fastq.gz")
				r2_file = os.path.join(raw_dir, f"{srr}_2.fastq.gz")
				
				if not (os.path.exists(r1_file) and os.path.exists(r2_file)):
					logger.error(f"Raw files not found for {srr}")
					update_stats('failed_processing')
				else:
					# Run FastQC
					success = True
					for fastq_file in [r1_file, r2_file]:
						fastqc_cmd = [args.fastqc_path, fastq_file, "-o", fastqc_dir, "--threads", "4", "--quiet"]
						fastqc_result = subprocess.run(fastqc_cmd, capture_output=True, text=True)
						if fastqc_result.returncode != 0:
							logger.warning(f"FastQC failed for {os.path.basename(fastq_file)}: {fastqc_result.stderr}")
							success = False
					
					if success:
						logger.info(f"FastQC completed for {srr}")
						update_stats('fastqc_completed')
						bbduk_queue.put([species, srr])  # Move to next stage
					else:
						update_stats('failed_processing')
						
			with worker_lock:
				active_fastqc_workers -= 1
				
			fastqc_queue.task_done()
			
		except queue.Empty:
			continue
		except Exception as e:
			logger.error(f"Error in FastQC worker: {str(e)}")
			with worker_lock:
				active_fastqc_workers -= 1

def bbduk_worker():
	"""Worker for BBduk cleaning"""
	global active_bbduk_workers
	logger = logging.getLogger(__name__)
	
	while not pipeline_shutdown.is_set():
		try:
			species, srr = bbduk_queue.get(timeout=1)
			if srr is None:  # Shutdown signal
				break
				
			with worker_lock:
				active_bbduk_workers += 1
				
			if check_cleaning_completed(srr):
				logger.info(f"BBduk cleaning already completed for {srr}, skipping")
				update_stats('cleaned')
				salmon_queue.put([species, srr])  # Move to next stage
			else:
				logger.info(f"Starting BBduk cleaning for {srr}")
				
				r1_file = os.path.join(raw_dir, f"{srr}_1.fastq.gz")
				r2_file = os.path.join(raw_dir, f"{srr}_2.fastq.gz")
				
				if not (os.path.exists(r1_file) and os.path.exists(r2_file)):
					logger.error(f"Raw files not found for {srr}")
					update_stats('failed_processing')
				else:
					clean_r1 = os.path.join(clean_dir, f"{srr}_clean_1.fastq.gz")
					clean_r2 = os.path.join(clean_dir, f"{srr}_clean_2.fastq.gz")
					
					for r_file, clean_r in [(r1_file, clean_r1), (r2_file, clean_r2)]:
						bbduk_cmd = [
							args.bbduk_path, f"in={r_file}", f"out={clean_r}", "threads=4"
						]

						bbduk_result = subprocess.run(bbduk_cmd, capture_output=True, text=True)
						if bbduk_result.returncode != 0:
							logger.error(f"Error cleaning {srr}: {bbduk_result.stderr}")
							update_stats('failed_processing')
						else:
							# Extract BBduk statistics from stderr
							reads_in, reads_out = 0, 0
							for line in bbduk_result.stderr.split('\n'):
								if 'Input:' in line and 'reads' in line:
									parts = line.split()
									for i, part in enumerate(parts):
										if part.isdigit() and 'reads' in parts[i+1]:
											reads_in = int(part)
											break
								elif 'Result:' in line and 'reads' in line:
									parts = line.split()
									for i, part in enumerate(parts):
										if part.isdigit() and 'reads' in parts[i+1]:
											reads_out = int(part)
											break
										
							retention_rate = (reads_out / reads_in * 100) if reads_in > 0 else 0
							logger.info(f"Successfully cleaned {srr} - Input: {reads_in:,} reads, Output: {reads_out:,} reads ({retention_rate:.1f}% retained)")

							# Remove original compressed files after successful cleaning
							try:
								if os.path.exists(r_file):
									os.remove(r_file)
								logger.debug(f"Removed original files for {srr}")
							except Exception as e:
								logger.warning(f"Could not remove original files for {srr}: {e}")

							update_stats('cleaned')
							salmon_queue.put([species, srr])  # Move to next stage
						
			with worker_lock:
				active_bbduk_workers -= 1
				
			bbduk_queue.task_done()
			
		except queue.Empty:
			continue
		except Exception as e:
			logger.error(f"Error in BBduk worker: {str(e)}")
			with worker_lock:
				active_bbduk_workers -= 1

def salmon_worker():
	"""Worker for Salmon quantification"""
	global active_salmon_workers
	logger = logging.getLogger(__name__)
	
	while not pipeline_shutdown.is_set():
		try:
			species, srr = salmon_queue.get(timeout=1)
			if srr is None:  # Shutdown signal
				break
				
			with worker_lock:
				active_salmon_workers += 1
				
			if check_quantification_completed(srr):
				logger.info(f"Salmon quantification already completed for {srr}, skipping")
				update_stats('quantified')
			else:
				# Wait for index to be built
				salmon_index_built.wait()
				
				logger.info(f"Starting Salmon quantification for {srr}")
				
				clean_r1 = os.path.join(clean_dir, f"{srr}_clean_1.fastq.gz")
				clean_r2 = os.path.join(clean_dir, f"{srr}_clean_2.fastq.gz")
				
				if not (os.path.exists(clean_r1) and os.path.exists(clean_r2)):
					logger.error(f"Clean files not found for {srr}")
					update_stats('failed_processing')
				else:
					sample_output = os.path.join(salmon_results_dir, srr)
					
					salmon_quant_cmd = [
						args.salmon_path, "quant", "-i", f"{salmon_index_dir}/{species}", "-l", "A",
						"-1", clean_r1, "-2", clean_r2,
						"-p", "4", "--validateMappings", "--seqBias", "--gcBias",
						"-o", sample_output
					]
					
					quant_result = subprocess.run(salmon_quant_cmd, capture_output=True, text=True)
					if quant_result.returncode != 0:
						logger.error(f"Error quantifying {srr}: {quant_result.stderr}")
						update_stats('failed_processing')
					else:
						# Extract Salmon statistics from stderr
						lib_type, mapping_rate = "Unknown", 0.0
						for line in quant_result.stderr.split('\n'):
							if 'library type' in line.lower() or 'detected library type' in line.lower():
								parts = line.split()
								for i, part in enumerate(parts):
									if part.upper() in ['ISF', 'ISR', 'IU', 'MSF', 'MSR', 'MU', 'OSF', 'OSR', 'OU', 'SF', 'SR', 'U']:
										lib_type = part.upper()
										break
							elif 'mapping rate' in line.lower():
								import re
								match = re.search(r'(\d+\.?\d*)%', line)
								if match:
									mapping_rate = float(match.group(1))
						
						logger.info(f"Successfully quantified {srr} - Library: {lib_type}, Mapping rate: {mapping_rate:.1f}%")
						update_stats('quantified')
						
						# Clean up cleaned files after successful quantification
						try:
							if os.path.exists(clean_r1):
								os.remove(clean_r1)
							if os.path.exists(clean_r2):
								os.remove(clean_r2)
							logger.debug(f"Removed cleaned files for {srr}")
						except Exception as e:
							logger.warning(f"Could not remove cleaned files for {srr}: {e}")
							
			with worker_lock:
				active_salmon_workers -= 1
				
			salmon_queue.task_done()
			
		except queue.Empty:
			continue
		except Exception as e:
			logger.error(f"Error in Salmon worker: {str(e)}")
			with worker_lock:
				active_salmon_workers -= 1

def get_available_slots():
	"""Calculate how many new workers can be started"""
	with worker_lock:
		total_active = active_fastqc_workers + active_bbduk_workers + active_salmon_workers
		return max(0, int(args.parallel_jobs) - total_active)

def build_salmon_index(species):
	"""Build Salmon index in a separate thread"""
	logger = logging.getLogger(__name__)
	
	# Check if index already exists
	index_info_file = os.path.join(salmon_index_dir, species, "info.json")
	if os.path.exists(index_info_file):
		logger.info("Salmon index already exists, skipping build")
		salmon_index_built.set()  # Signal that index is ready
		return
	
	logger.info(f"Building Salmon index for species: {species}")
	
	if not os.path.exists(args.transcripts):
		logger.error(f"Transcriptome file not found: {args.transcripts}")
		sys.exit(1)
	
	# Build index command
	salmon_index_cmd = [args.salmon_path, "index", "-t", args.transcripts, "-i", f"{salmon_index_dir}/{species}"]
	
	# Add decoy if provided
	if args.decoy:
		if not os.path.exists(args.decoy):
			logger.error(f"Decoy file not found: {args.decoy}")
			sys.exit(1)
		
		# Create combined transcriptome + decoy file
		combined_file = os.path.join(salmon_index_dir, f"{species}_combined.fa")
		logger.info(f"Creating combined transcriptome + decoy file: {combined_file}")
		
		with open(combined_file, 'w') as out_f:
			# Copy transcriptome
			with open(args.transcripts, 'r') as trans_f:
				out_f.write(trans_f.read())
			# Copy decoy
			with open(args.decoy, 'r') as decoy_f:
				out_f.write(decoy_f.read())
		
		# Create decoy names file
		decoy_names_file = os.path.join(salmon_index_dir, f"{species}_decoys.txt")
		with open(decoy_names_file, 'w') as decoy_f:
			with open(args.decoy, 'r') as input_f:
				for line in input_f:
					if line.startswith('>'):
						decoy_name = line[1:].split()[0]  # Remove '>' and get first part
						decoy_f.write(f"{decoy_name}\n")
		
		# Update command to use combined file and decoy names
		salmon_index_cmd = [
			args.salmon_path, "index", "-t", combined_file, "-i", f"{salmon_index_dir}/{species}",
			"-d", decoy_names_file
		]
		logger.info(f"Using decoy sequences from: {args.decoy}")
	
	index_result = subprocess.run(salmon_index_cmd, capture_output=True, text=True)
	if index_result.returncode != 0:
		logger.error(f"Error building Salmon index: {index_result.stderr}")
		sys.exit(1)
	
	logger.info("Salmon index successfully created")
	salmon_index_built.set()  # Signal that index is ready

# Main pipeline execution
logger.info("=" * 80)
logger.info("ORDERED PIPELINE RNA-SEQ PROCESSING")
logger.info("=" * 80)
logger.info(f"Total species: {len(sra_ids)}")
logger.info(f"Parallel jobs: {args.parallel_jobs}")
logger.info(f"Base URL: {args.base_url}")
logger.info("Pipeline structure:")
logger.info("  Download: Sequential (one sample at a time)")
logger.info("  Processing: Parallel per sample (FastQC → BBduk → Salmon)")
logger.info("  Threads per tool: FastQC=4, BBduk=4, Salmon=4")
logger.info(f"  Total concurrent samples: Up to {args.parallel_jobs}")

# Check for existing results and show resume information
logger.info("CHECKING EXISTING RESULTS")

# Start building Salmon index in background
index_thread = threading.Thread(target=build_salmon_index, args=(list(sra_ids.keys())[0],))
index_thread.daemon = True
index_thread.start()

# Start workers for parallel processing
logger.info("Starting worker threads")
logger.info(f"Architecture: {args.parallel_jobs} parallel samples × 3 stages (FastQC/BBduk/Salmon) × 4 threads each")
workers = []

# Start worker threads for each stage
for _ in range(int(args.parallel_jobs)):
	# FastQC workers
	worker = threading.Thread(target=fastqc_worker)
	worker.daemon = True
	worker.start()
	workers.append(worker)
	
	# BBduk workers  
	worker = threading.Thread(target=bbduk_worker)
	worker.daemon = True
	worker.start()
	workers.append(worker)
	
	# Salmon workers
	worker = threading.Thread(target=salmon_worker)
	worker.daemon = True
	worker.start()
	workers.append(worker)

logger.info(f"Started {len(workers)} worker threads")

# Main pipeline execution - Sequential download with ordered parallel processing
logger.info("Starting sequential download and ordered parallel processing")

for i, (species, srr_list) in enumerate(sra_ids.items(), 1):
	for srr in srr_list:
		# Skip if already fully processed
		if check_quantification_completed(srr):
			logger.info(f"{srr} already fully processed, skipping")
			continue

		logger.info(f"Processing sample {i}/{len(sra_ids)}: {srr}")

		# Se já existe FastQC, não baixa, só manda para o BBduk (se passar critério)
		if check_fastqc_completed(srr):
			logger.info(f"{srr} FastQC already exists, skipping download and FastQC step")
			if not check_cleaning_completed(srr):
				if fastqc_passed_quality(srr):
					bbduk_queue.put([species, srr])
					logger.info(f"{srr} added to BBduk queue (FastQC passed quality)")
				else:
					logger.warning(f"{srr} did not pass FastQC quality threshold (>80% good sequences), skipping BBduk and Salmon")
			elif not check_quantification_completed(srr):
				if check_cleaning_completed(srr):
					salmon_queue.put([species, srr])
					logger.info(f"{srr} added to Salmon queue (BBduk completed)")
			continue

		# Download sequencialmente (apenas se não existe FastQC)
		if download_sample(srr):
			if not check_fastqc_completed(srr):
				fastqc_queue.put([species, srr])
				logger.info(f"{srr} added to FastQC queue")
			elif not check_cleaning_completed(srr):
				if fastqc_passed_quality(srr):
					bbduk_queue.put([species, srr])
					logger.info(f"{srr} added to BBduk queue (FastQC passed quality)")
				else:
					logger.warning(f"{srr} did not pass FastQC quality threshold (>80% good sequences), skipping BBduk and Salmon")
			elif not check_quantification_completed(srr):
				if check_cleaning_completed(srr):
					salmon_queue.put([species, srr])
					logger.info(f"{srr} added to Salmon queue (BBduk completed)")
		else:
			logger.error(f"{srr} download failed, skipping processing")

# Wait for all queues to be processed
logger.info("Waiting for all processing to complete")

# Monitor progress
while (not fastqc_queue.empty() or 
		not bbduk_queue.empty() or 
		not salmon_queue.empty() or
		active_fastqc_workers > 0 or 
		active_bbduk_workers > 0 or 
		active_salmon_workers > 0):
	
	with worker_lock:
		total_queued = fastqc_queue.qsize() + bbduk_queue.qsize() + salmon_queue.qsize()
		total_active = active_fastqc_workers + active_bbduk_workers + active_salmon_workers
		logger.info(f"Progress - Queued: {total_queued} samples | Active: {total_active} workers")
		logger.debug(f"Breakdown - FastQC: {fastqc_queue.qsize()}q/{active_fastqc_workers}w, BBduk: {bbduk_queue.qsize()}q/{active_bbduk_workers}w, Salmon: {salmon_queue.qsize()}q/{active_salmon_workers}w")
	
	time.sleep(5)  # Check every 5 seconds

logger.info("All processing completed")
	
# Shutdown workers
logger.info("Shutting down worker threads")
pipeline_shutdown.set()

# Send shutdown signals to queues
for _ in range(int(args.parallel_jobs)):
	fastqc_queue.put(None, None)
	bbduk_queue.put(None, None)
	salmon_queue.put(None, None)

# Wait for workers to finish
for worker in workers:
	worker.join(timeout=10)

# Wait for index building to complete
index_thread.join()

logger.info("MERGING QUANTIFICATION RESULTS")

# Merge quantification results
if processing_stats['quantified'] > 0:
	logger.info("Merging quantification results")
	
	# Find all directories with quant.sf files
	quant_dirs = []
	for item in os.listdir(salmon_results_dir):
		quant_file = os.path.join(salmon_results_dir, item, "quant.sf")
		if os.path.isfile(quant_file):
			quant_dirs.append(os.path.join(salmon_results_dir, item))
	
	if quant_dirs:
		# Merge TPM values
		tpm_merge_cmd = [
			args.salmon_path, "quantmerge",
			"--column", "TPM",
			"--output", os.path.join(outdir, "merged_TPM.tsv"),
			"--quants"
		] + quant_dirs
		
		tpm_result = subprocess.run(tpm_merge_cmd, capture_output=True, text=True)
		if tpm_result.returncode != 0:
			logger.error(f"Error merging TPM: {tpm_result.stderr}")
		else:
			logger.info("TPM merge completed successfully")
		
		# Merge count values
		counts_merge_cmd = [
			args.salmon_path, "quantmerge",
			"--column", "NumReads",
			"--output", os.path.join(outdir, "merged_counts.tsv"),
			"--quants"
		] + quant_dirs
		
		counts_result = subprocess.run(counts_merge_cmd, capture_output=True, text=True)
		if counts_result.returncode != 0:
			logger.error(f"Error merging counts: {counts_result.stderr}")
		else:
			logger.info("Counts merge completed successfully")

logger.info("FINAL CLEANUP")

# Clean up remaining raw and clean directories
logger.info("Removing remaining temporary files")
try:
	if os.path.exists(raw_dir):
		raw_size = sum(os.path.getsize(os.path.join(dirpath, filename))
						for dirpath, dirnames, filenames in os.walk(raw_dir)
						for filename in filenames) / (1024*1024)  # Size in MB
		
		shutil.rmtree(raw_dir)
		logger.info(f"Removed raw data directory ({raw_size:.1f} MB freed)")
		
except Exception as e:
	logger.warning(f"Error during cleanup: {e}")

logger.info("=" * 80)
logger.info("PIPELINE COMPLETED!")
logger.info("=" * 80)

# Final statistics
logger.info("Processing Statistics:")
logger.info(f"  Downloaded: {processing_stats['downloaded']}/{len(sra_ids)}")
logger.info(f"  FastQC completed: {processing_stats['fastqc_completed']}/{len(sra_ids)}")
logger.info(f"  Cleaned: {processing_stats['cleaned']}/{len(sra_ids)}")
logger.info(f"  Quantified: {processing_stats['quantified']}/{len(sra_ids)}")
logger.info(f"  Failed downloads: {processing_stats['failed_downloads']}")
logger.info(f"  Failed processing: {processing_stats['failed_processing']}")

logger.info(f"Results directory: {outdir}")
if processing_stats['downloaded'] > 0:
	logger.info(f"  - FastQC results: {fastqc_dir}")
	logger.info(f"  - Salmon quantification results: {salmon_results_dir}")
	if os.path.exists(os.path.join(outdir, "merged_TPM.tsv")):
		logger.info(f"  - Merged TPM results: {os.path.join(outdir, 'merged_TPM.tsv')}")
	if os.path.exists(os.path.join(outdir, "merged_counts.tsv")):
		logger.info(f"  - Merged counts results: {os.path.join(outdir, 'merged_counts.tsv')}")

# Exit with appropriate code
total_failures = processing_stats['failed_downloads'] + processing_stats['failed_processing']
sys.exit(0 if total_failures == 0 else 1)
