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
	python3 preproc_pipeline.py --sra_list samples.txt --transcripts transcriptome.fa --decoy decoy.fa --salmon_path /path/to/salmon --fastqc_path /path/to/fastqc --bbduk_path /path/to/bbduk --config config.yaml --outdir output --parallel_jobs 4
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
from contextlib import contextmanager

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--sra_list", required=True, help="File with the IDs (SRR)")
parser.add_argument("--transcripts", required=True, help="Path to transcriptome FASTA file (applied to all species), or TSV file with 'species\\tpath' per line for per-species transcriptomes")
parser.add_argument("--config", default="config.yaml", help="Path to configuration YAML file")
parser.add_argument("--outdir", default="pipeline_output", help="Output directory")
parser.add_argument("--salmon_path", default="salmon", help="Path to Salmon")
parser.add_argument("--fastqc_path", default="fastqc", help="Path to FastQC")
parser.add_argument("--bbduk_path", default="bbduk.sh", help="Path to BBduk")
parser.add_argument("--parallel_jobs", default="5", help="Maximum number of parallel jobs")
parser.add_argument("--cores", default="4", help="Number of threads per tool (FastQC, BBduk, Salmon)")
parser.add_argument("--decoy", help="Path to decoy sequences FASTA file (optional, applied to all species), or TSV file with 'species\\tpath' per line for per-species decoys")
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
error_dir = os.path.join(outdir, "error")
os.makedirs(error_dir, exist_ok=True)

# Initialize logging
logger = setup_logging(outdir)
logger.info("Pipeline initialized successfully")
logger.info(f"Output directory: {outdir}")

# Read IDs and fill list with SRRs
sra_ids = {}
with open(args.sra_list) as f:
	for line in f:
		parts = line.strip().split()
		if len(parts) < 2:
			continue
		species = parts[0]
		srr = parts[-1]  # last column is always the SRR ID (handles 2 or 3+ column formats)
		if species not in sra_ids:
			sra_ids[species] = [srr]
		else:
			sra_ids[species].append(srr)

def _load_tsv_map(filepath):
	"""Load a TSV file with 'species\\tpath' lines into a dict."""
	result = {}
	with open(filepath, 'r') as _f:
		for _line in _f:
			_line = _line.strip()
			if _line and not _line.startswith('#'):
				_parts = _line.split('\t', 1)
				if len(_parts) == 2:
					result[_parts[0]] = _parts[1]
	return result

def _is_tsv_map(filepath):
	"""Return True if file looks like a species->path TSV (not a FASTA)."""
	with open(filepath, 'r') as _f:
		_first = _f.readline().strip()
	return '\t' in _first and not _first.startswith('>')

# Build transcripts map: {species: path}
transcripts_map = {}
if _is_tsv_map(args.transcripts):
	multi_species_mode = True
	transcripts_map = _load_tsv_map(args.transcripts)
	for _sp in sra_ids:
		if _sp not in transcripts_map:
			print(f"No transcript file found for species '{_sp}' in transcripts mapping file")
			sys.exit(1)
else:
	multi_species_mode = False
	for _sp in sra_ids:
		transcripts_map[_sp] = args.transcripts

# Build decoys map: {species: path | None}
decoys_map = {}
if args.decoy and _is_tsv_map(args.decoy):
	_decoy_tsv = _load_tsv_map(args.decoy)
	for _sp in sra_ids:
		decoys_map[_sp] = _decoy_tsv.get(_sp)  # None if not present
elif args.decoy:
	for _sp in sra_ids:
		decoys_map[_sp] = args.decoy
else:
	for _sp in sra_ids:
		decoys_map[_sp] = None

# Create per-species salmon results directories when in multi-species mode
species_salmon_results_dir = {}
if multi_species_mode:
	for _sp in sra_ids:
		_sp_dir = os.path.join(outdir, _sp, "salmon_results")
		os.makedirs(_sp_dir, exist_ok=True)
		species_salmon_results_dir[_sp] = _sp_dir
else:
	for _sp in sra_ids:
		species_salmon_results_dir[_sp] = salmon_results_dir

# Global variables for tracking progress
sample_queue = queue.Queue()
fastqc_queue = queue.Queue()
bbduk_queue = queue.Queue()
salmon_queue = queue.Queue()
salmon_index_built = {species: threading.Event() for species in sra_ids}
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

# Global core pool: limits total threads used across all concurrent tool processes
class CorePool:
    """Thread-safe pool that limits total cores used across all running tools."""
    def __init__(self, total_cores):
        self._available = total_cores
        self._total = total_cores
        self._cond = threading.Condition()

    def acquire(self, n):
        with self._cond:
            while self._available < n:
                self._cond.wait()
            self._available -= n

    def release(self, n):
        with self._cond:
            self._available += n
            self._cond.notify_all()

    @contextmanager
    def use(self, n):
        self.acquire(n)
        try:
            yield n
        finally:
            self.release(n)

total_cores = int(args.cores)
parallel_jobs = int(args.parallel_jobs)
cores_per_job = max(1, total_cores // parallel_jobs)
core_pool = CorePool(total_cores)

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
			# Example: <tr><td>Total Sequences</td><td>23400214</td></tr>
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

def check_quantification_completed(srr, species):
	"""Check if Salmon quantification is already completed for a sample"""
	sp_results_dir = species_salmon_results_dir[species]
	quant_log = os.path.join(sp_results_dir, srr, "logs", "salmon_quant.log")
	quant_file = os.path.join(sp_results_dir, srr, "quant.sf")

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
					# Save error
					error_file = os.path.join(error_dir, f"{srr}_fastqc_error.txt")
					try:
						with open(error_file, 'w') as ef:
							ef.write("Raw files not found for FastQC\n")
					except Exception as e:
						logger.warning(f"Could not write FastQC error file for {srr}: {e}")
					# Remove raw data
					for f in [r1_file, r2_file]:
						try:
							if os.path.exists(f):
								os.remove(f)
								logger.debug(f"Removed raw file after FastQC error: {f}")
						except Exception as e:
							logger.warning(f"Could not remove raw file {f} after FastQC error: {e}")
				else:
					# Run FastQC
					success = True
					fastqc_errors = ""
					for fastq_file in [r1_file, r2_file]:
						fastqc_cmd = [args.fastqc_path, fastq_file, "-o", fastqc_dir, "--threads", str(cores_per_job), "--quiet"]
						with core_pool.use(cores_per_job):
							fastqc_result = subprocess.run(fastqc_cmd, capture_output=True, text=True)
						if fastqc_result.returncode != 0:
							logger.warning(f"FastQC failed for {os.path.basename(fastq_file)}: {fastqc_result.stderr}")
							success = False
							fastqc_errors += f"FastQC failed for {os.path.basename(fastq_file)}: {fastqc_result.stderr}\n"
					if success:
						logger.info(f"FastQC completed for {srr}")
						update_stats('fastqc_completed')
						bbduk_queue.put([species, srr])  # Move to next stage
					else:
						update_stats('failed_processing')
						# Save error
						error_file = os.path.join(error_dir, f"{srr}_fastqc_error.txt")
						try:
							with open(error_file, 'w') as ef:
								ef.write(fastqc_errors)
						except Exception as e:
							logger.warning(f"Could not write FastQC error file for {srr}: {e}")
						# Remove raw data
						for f in [r1_file, r2_file]:
							try:
								if os.path.exists(f):
									os.remove(f)
									logger.debug(f"Removed raw file after FastQC error: {f}")
							except Exception as e:
								logger.warning(f"Could not remove raw file {f} after FastQC error: {e}")
						
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
					# Save error
					error_file = os.path.join(error_dir, f"{srr}_bbduk_error.txt")
					try:
						with open(error_file, 'w') as ef:
							ef.write("Raw files not found for BBduk\n")
					except Exception as e:
						logger.warning(f"Could not write BBduk error file for {srr}: {e}")
					# Remove raw data
					for f in [r1_file, r2_file]:
						try:
							if os.path.exists(f):
								os.remove(f)
								logger.debug(f"Removed raw file after BBduk error: {f}")
						except Exception as e:
							logger.warning(f"Could not remove raw file {f} after BBduk error: {e}")
				else:
					clean_r1 = os.path.join(clean_dir, f"{srr}_clean_1.fastq.gz")
					clean_r2 = os.path.join(clean_dir, f"{srr}_clean_2.fastq.gz")

					success_r1 = False
					success_r2 = False
					for idx, (r_file, clean_r) in enumerate([(r1_file, clean_r1), (r2_file, clean_r2)]):
						bbduk_cmd = [
							args.bbduk_path, f"in={r_file}", f"out={clean_r}", f"threads={cores_per_job}"
						]

						with core_pool.use(cores_per_job):
							bbduk_result = subprocess.run(bbduk_cmd, capture_output=True, text=True)
						if bbduk_result.returncode != 0:
							logger.error(f"Error cleaning {srr} ({'R1' if idx==0 else 'R2'}): {bbduk_result.stderr}")
							update_stats('failed_processing')
							# Save error
							error_file = os.path.join(error_dir, f"{srr}_bbduk_error.txt")
							try:
								with open(error_file, 'a') as ef:
									ef.write(f"Error cleaning {srr} ({'R1' if idx==0 else 'R2'}): {bbduk_result.stderr}\n")
							except Exception as e:
								logger.warning(f"Could not write BBduk error file for {srr}: {e}")
							# Remove raw data
							try:
								if os.path.exists(r_file):
									os.remove(r_file)
									logger.debug(f"Removed raw file after BBduk error: {r_file}")
							except Exception as e:
								logger.warning(f"Could not remove raw file {r_file} after BBduk error: {e}")
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
							logger.info(f"Successfully cleaned {srr} ({'R1' if idx==0 else 'R2'}) - Input: {reads_in:,} reads, Output: {reads_out:,} reads ({retention_rate:.1f}% retained)")

							# Remove original compressed files after successful cleaning
							try:
								if os.path.exists(r_file):
									os.remove(r_file)
								logger.debug(f"Removed original files for {srr} ({'R1' if idx==0 else 'R2'})")
							except Exception as e:
								logger.warning(f"Could not remove original files for {srr} ({'R1' if idx==0 else 'R2'}): {e}")

							if idx == 0:
								success_r1 = True
							else:
								success_r2 = True

					# Only updates stats and adds to salmon_queue if both were cleaned successfully
					if success_r1 and success_r2:
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
				
			if check_quantification_completed(srr, species):
				logger.info(f"Salmon quantification already completed for {srr}, skipping")
				update_stats('quantified')
				# Remove raw data
				r1_file = os.path.join(raw_dir, f"{srr}_1.fastq.gz")
				r2_file = os.path.join(raw_dir, f"{srr}_2.fastq.gz")
				for f in [r1_file, r2_file]:
					try:
						if os.path.exists(f):
							os.remove(f)
							logger.debug(f"Removed raw file after Salmon: {f}")
					except Exception as e:
						logger.warning(f"Could not remove raw file {f} after Salmon: {e}")
			else:
				# Wait for index to be built for this species
				salmon_index_built[species].wait()

				logger.info(f"Starting Salmon quantification for {srr}")

				clean_r1 = os.path.join(clean_dir, f"{srr}_clean_1.fastq.gz")
				clean_r2 = os.path.join(clean_dir, f"{srr}_clean_2.fastq.gz")

				r1_file = os.path.join(raw_dir, f"{srr}_1.fastq.gz")
				r2_file = os.path.join(raw_dir, f"{srr}_2.fastq.gz")

				if not (os.path.exists(clean_r1) and os.path.exists(clean_r2)):
					logger.error(f"Clean files not found for {srr}")
					update_stats('failed_processing')
					# Remove raw data
					for f in [r1_file, r2_file]:
						try:
							if os.path.exists(f):
								os.remove(f)
								logger.debug(f"Removed raw file after Salmon error: {f}")
						except Exception as e:
							logger.warning(f"Could not remove raw file {f} after Salmon error: {e}")
				else:
					sample_output = os.path.join(species_salmon_results_dir[species], srr)

					salmon_quant_cmd = [
						args.salmon_path, "quant", "-i", f"{salmon_index_dir}/{species}", "-l", "A",
						"-1", clean_r1, "-2", clean_r2,
						"-p", str(cores_per_job), "--validateMappings", "--seqBias", "--gcBias", "--dumpEq", "--numGibbsSamples", "200",
						"-o", sample_output
					]

					with core_pool.use(cores_per_job):
						quant_result = subprocess.run(salmon_quant_cmd, capture_output=True, text=True)
					if quant_result.returncode != 0:
						logger.error(f"Error quantifying {srr}: {quant_result.stderr}")
						update_stats('failed_processing')
						# Write error to specific file
						error_file = os.path.join(error_dir, f"{srr}_salmon_error.txt")
						try:
							with open(error_file, 'w') as ef:
								ef.write(quant_result.stderr)
						except Exception as e:
							logger.warning(f"Could not write error file for {srr}: {e}")
						# Remove raw data
						for f in [r1_file, r2_file]:
							try:
								if os.path.exists(f):
									os.remove(f)
									logger.debug(f"Removed raw file after Salmon error: {f}")
							except Exception as e:
								logger.warning(f"Could not remove raw file {f} after Salmon error: {e}")
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
						# Remove raw data
						for f in [r1_file, r2_file]:
							try:
								if os.path.exists(f):
									os.remove(f)
									logger.debug(f"Removed raw file after Salmon: {f}")
							except Exception as e:
								logger.warning(f"Could not remove raw file {f} after Salmon: {e}")
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

def build_salmon_index(species, transcript_path):
	"""Build Salmon index in a separate thread"""
	logger = logging.getLogger(__name__)

	# Check if index already exists
	index_info_file = os.path.join(salmon_index_dir, species, "info.json")
	if os.path.exists(index_info_file):
		logger.info(f"Salmon index for '{species}' already exists, skipping build")
		salmon_index_built[species].set()  # Signal that index is ready
		return

	logger.info(f"Building Salmon index for species: {species}")

	if not os.path.exists(transcript_path):
		logger.error(f"Transcriptome file not found for species '{species}': {transcript_path}")
		sys.exit(1)

	# Build index command
	salmon_index_cmd = [args.salmon_path, "index", "-t", transcript_path, "-i", f"{salmon_index_dir}/{species}"]

	# Add decoy if provided for this species
	decoy_path = decoys_map.get(species)
	if decoy_path:
		if not os.path.exists(decoy_path):
			logger.error(f"Decoy file not found for species '{species}': {decoy_path}")
			sys.exit(1)

		# Create combined transcriptome + decoy file
		combined_file = os.path.join(salmon_index_dir, f"{species}_combined.fa")
		logger.info(f"Creating combined transcriptome + decoy file: {combined_file}")

		with open(combined_file, 'w') as out_f:
			# Copy transcriptome
			with open(transcript_path, 'r') as trans_f:
				out_f.write(trans_f.read())
			# Copy decoy
			with open(decoy_path, 'r') as decoy_f:
				out_f.write(decoy_f.read())

		# Create decoy names file
		decoy_names_file = os.path.join(salmon_index_dir, f"{species}_decoys.txt")
		with open(decoy_names_file, 'w') as decoy_f:
			with open(decoy_path, 'r') as input_f:
				for line in input_f:
					if line.startswith('>'):
						decoy_name = line[1:].split()[0]  # Remove '>' and get first part
						decoy_f.write(f"{decoy_name}\n")

		# Update command to use combined file and decoy names
		salmon_index_cmd = [
			args.salmon_path, "index", "-t", combined_file, "-i", f"{salmon_index_dir}/{species}",
			"-d", decoy_names_file
		]
		logger.info(f"Using decoy sequences from: {decoy_path}")

	index_result = subprocess.run(salmon_index_cmd, capture_output=True, text=True)
	if index_result.returncode != 0:
		logger.error(f"Error building Salmon index for '{species}': {index_result.stderr}")
		sys.exit(1)

	logger.info(f"Salmon index for '{species}' successfully created")
	salmon_index_built[species].set()  # Signal that index is ready

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
logger.info("  Core pool: {total_cores} total cores, {cores_per_job} cores/job ({parallel_jobs} parallel jobs)")
logger.info(f"  Total concurrent samples: Up to {args.parallel_jobs}")

# Check for existing results and show resume information
logger.info("CHECKING EXISTING RESULTS")

# Start building Salmon index for each species in background
logger.info(f"Starting Salmon index building for {len(transcripts_map)} species")
index_threads = []
for _species, _transcript_path in transcripts_map.items():
	_t = threading.Thread(target=build_salmon_index, args=(_species, _transcript_path))
	_t.daemon = True
	_t.start()
	index_threads.append(_t)

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
		if check_quantification_completed(srr, species):
			logger.info(f"{srr} already fully processed, skipping")
			continue

		logger.info(f"Processing sample {i}/{len(sra_ids)}: {srr}")

		# If FastQC already exists, skip download and send directly to BBduk (if criteria are met)
		if check_fastqc_completed(srr):
			logger.info(f"{srr} FastQC already exists, skipping download and FastQC step")
			if not check_cleaning_completed(srr):
				if fastqc_passed_quality(srr):
					bbduk_queue.put([species, srr])
					logger.info(f"{srr} added to BBduk queue (FastQC passed quality)")
				else:
					logger.warning(f"{srr} did not pass FastQC quality threshold (>80% good sequences), skipping BBduk and Salmon")
			elif not check_quantification_completed(srr, species):
				if check_cleaning_completed(srr):
					salmon_queue.put([species, srr])
					logger.info(f"{srr} added to Salmon queue (BBduk completed)")
			continue

		# Download sequentially (only if FastQC does not exist)
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
			elif not check_quantification_completed(srr, species):
				if check_cleaning_completed(srr):
					salmon_queue.put([species, srr])
					logger.info(f"{srr} added to Salmon queue (BBduk completed)")
		else:
			logger.error(f"{srr} download failed, skipping processing")

# Wait for all queues to be processed
logger.info("Waiting for all processing to complete")

# Monitor progress
prev_queued = None
prev_active = None
while (not fastqc_queue.empty() or 
		not bbduk_queue.empty() or 
		not salmon_queue.empty() or
		active_fastqc_workers > 0 or 
		active_bbduk_workers > 0 or 
		active_salmon_workers > 0):

	with worker_lock:
		total_queued = fastqc_queue.qsize() + bbduk_queue.qsize() + salmon_queue.qsize()
		total_active = active_fastqc_workers + active_bbduk_workers + active_salmon_workers
		if total_queued != prev_queued or total_active != prev_active:
			logger.info(f"Progress - Queued: {total_queued} samples | Active: {total_active} workers")
			logger.debug(f"Breakdown - FastQC: {fastqc_queue.qsize()}q/{active_fastqc_workers}w, BBduk: {bbduk_queue.qsize()}q/{active_bbduk_workers}w, Salmon: {salmon_queue.qsize()}q/{active_salmon_workers}w")
			prev_queued = total_queued
			prev_active = total_active
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

# Wait for all index building threads to complete
for _t in index_threads:
	_t.join()

logger.info("MERGING QUANTIFICATION RESULTS")

# Merge quantification results per species
if processing_stats['quantified'] > 0:
	logger.info("Merging quantification results")

	for _sp, _sp_results_dir in species_salmon_results_dir.items():
		_merge_outdir = os.path.join(outdir, _sp) if multi_species_mode else outdir

		# Find all directories with quant.sf files for this species
		quant_dirs = []
		if os.path.isdir(_sp_results_dir):
			for item in os.listdir(_sp_results_dir):
				quant_file = os.path.join(_sp_results_dir, item, "quant.sf")
				if os.path.isfile(quant_file):
					quant_dirs.append(os.path.join(_sp_results_dir, item))

		if not quant_dirs:
			continue

		_sp_label = f" [{_sp}]" if multi_species_mode else ""

		# Merge TPM values
		tpm_merge_cmd = [
			args.salmon_path, "quantmerge",
			"--column", "TPM",
			"--output", os.path.join(_merge_outdir, "merged_TPM.tsv"),
			"--quants"
		] + quant_dirs

		tpm_result = subprocess.run(tpm_merge_cmd, capture_output=True, text=True)
		if tpm_result.returncode != 0:
			logger.error(f"Error merging TPM{_sp_label}: {tpm_result.stderr}")
		else:
			logger.info(f"TPM merge completed successfully{_sp_label}")

		# Merge count values
		counts_merge_cmd = [
			args.salmon_path, "quantmerge",
			"--column", "NumReads",
			"--output", os.path.join(_merge_outdir, "merged_counts.tsv"),
			"--quants"
		] + quant_dirs

		counts_result = subprocess.run(counts_merge_cmd, capture_output=True, text=True)
		if counts_result.returncode != 0:
			logger.error(f"Error merging counts{_sp_label}: {counts_result.stderr}")
		else:
			logger.info(f"Counts merge completed successfully{_sp_label}")

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
	if multi_species_mode:
		for _sp, _sp_dir in species_salmon_results_dir.items():
			logger.info(f"  - Salmon results [{_sp}]: {_sp_dir}")
			_sp_outdir = os.path.join(outdir, _sp)
			if os.path.exists(os.path.join(_sp_outdir, "merged_TPM.tsv")):
				logger.info(f"  - Merged TPM [{_sp}]: {os.path.join(_sp_outdir, 'merged_TPM.tsv')}")
			if os.path.exists(os.path.join(_sp_outdir, "merged_counts.tsv")):
				logger.info(f"  - Merged counts [{_sp}]: {os.path.join(_sp_outdir, 'merged_counts.tsv')}")
	else:
		logger.info(f"  - Salmon quantification results: {salmon_results_dir}")
		if os.path.exists(os.path.join(outdir, "merged_TPM.tsv")):
			logger.info(f"  - Merged TPM results: {os.path.join(outdir, 'merged_TPM.tsv')}")
		if os.path.exists(os.path.join(outdir, "merged_counts.tsv")):
			logger.info(f"  - Merged counts results: {os.path.join(outdir, 'merged_counts.tsv')}")

# Exit with appropriate code
total_failures = processing_stats['failed_downloads'] + processing_stats['failed_processing']
sys.exit(0 if total_failures == 0 else 1)
