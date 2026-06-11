#!/usr/bin/env python3
import os
import sys
import yaml
import logging
import time
import subprocess

# Import our custom step modules
import step1_trimming
import step2_fastq_to_bam
import step3_umi_align_consensus
import step4_base_recaliberation
import step5_variant_calling

def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: python pipeline_runner.py <working_directory> <sample_info> <config.yaml>")

    # Make paths absolute right away to avoid working directory relative issues
    project_dir = os.path.abspath(sys.argv[1])
    sample_info_file = os.path.abspath(sys.argv[2])
    config_file = os.path.abspath(sys.argv[3])

    # Setup the structured pipeline logging directory and name
    log_name = f"pipeline_{int(time.time())}.log"
    init_logging(project_dir, log_name)
    logger = logging.getLogger('pipeline')
    
    logger.info("==================================================")
    logger.info("Initializing ctDNA Pipeline")
    logger.info(f"Working Dir: {project_dir}")
    logger.info("==================================================")

    # Core configurations
    config = read_config(config_file)
    sample_ids = read_sample_info(sample_info_file)

    # Setup required working cluster infrastructure directories
    for folder in ["job_files", "job_output", "job_error", "bam", "trimmed_fastq"]:
        os.makedirs(os.path.join(project_dir, folder), exist_ok=True)

    # -------------------------------------------------------------------------
    # STEP 1: READ TRIMMING
    # -------------------------------------------------------------------------
    if config.get("run_steps", {}).get("read_trimming", False):
        logger.info("--- Starting Step 1: Read Trimming ---")
        step1_jobs = step1_trimming.run_step(project_dir, sample_info_file, config_file, sample_ids)
        wait_for_jobs(step1_jobs, logger)
        logger.info("Step 1: Read Trimming completed successfully.")
    else:
        logger.info("--- Skipping Step 1: Read Trimming (Disabled in config) ---")

    # -------------------------------------------------------------------------
    # STEP 2: FASTQ TO BAM ALIGNMENT
    # -------------------------------------------------------------------------
    if config.get("run_steps", {}).get("merge_fastq", False) or config.get("run_steps", {}).get("fastq_to_bam", False):
        logger.info("--- Starting Step 2: Fastq To BAM & ZipperBams ---")
        step2_jobs = step2_fastq_to_bam.run_step(project_dir, sample_info_file, config_file)
        wait_for_jobs(step2_jobs, logger)
        logger.info("Step 2: Fastq to BAM alignment completed successfully.")
    else:
        logger.info("--- Skipping Step 2: Fastq To BAM & ZipperBams (Disabled in config) ---")

    # -------------------------------------------------------------------------
    # STEP 3: UMI ALIGNMENT & MOLECULAR CONSENSUS
    # -------------------------------------------------------------------------
    if config.get("run_steps", {}).get("align_consensus_umis", False):
        logger.info("--- Starting Step 3: UMI Grouping and Consensus Calling ---")
        step3_jobs = step3_umi_align_consensus.run_step(project_dir, sample_info_file, config_file, sample_ids)
        wait_for_jobs(step3_jobs, logger)
        logger.info("Step 3: UMI Tracking & Consensus completed successfully.")
    else:
        logger.info("--- Skipping Step 3: UMI Alignment & Consensus (Disabled in config) ---")

    # -------------------------------------------------------------------------
    # STEP 4: BASE QUALITY SCORE RECALIBRATION (BQSR)
    # -------------------------------------------------------------------------
    if config.get("run_steps", {}).get("base_recalibration", False):
        logger.info("--- Starting Step 4: Base Quality Score Recalibration ---")
        step4_jobs = step4_base_recaliberation.run_step(project_dir, sample_info_file, config_file, sample_ids)
        wait_for_jobs(step4_jobs, logger)
        logger.info("Step 4: BQSR processing completed successfully.")
    else:
        logger.info("--- Skipping Step 4: Base Quality Score Recalibration (Disabled in config) ---")

    # -------------------------------------------------------------------------
    # STEP 5: VARIANT CALLING & ANNOTATION (MUTECT2 / FREEBAYES / VEP / FUNCOTATOR)
    # -------------------------------------------------------------------------
    if config.get("run_steps", {}).get("variant_and_copy_number", False) or config.get("run_steps", {}).get("variant_calling", False):
        logger.info("--- Starting Step 5: Variant Calling, Normalization & Annotation ---")
        step5_jobs = step5_variant_calling.run_step(project_dir, sample_info_file, config_file, sample_ids)
        wait_for_jobs(step5_jobs, logger)
        logger.info("Step 5: Variant Calling & Annotation completed successfully.")
    else:
        logger.info("--- Skipping Step 5: Variant Calling & Annotation (Disabled in config) ---")
        
    logger.info("==================================================")
    logger.info("Pipeline Complete.")
    logger.info("==================================================")


def wait_for_jobs(job_ids, logger):
    """
    Tracks cluster progress specifically matching SLURM Job IDs to names.
    Logs explicitly upon initialization, individual job completion, or failures.
    """
    if not job_ids:
        logger.warning("No jobs submitted or running for this phase.")
        return

    job_ids = [j for j in job_ids if j]
    
    # Pre-populate tracking map using scontrol to verify submission names immediately
    tracked_jobs = {}
    for j_id in job_ids:
        try:
            res = subprocess.run(["scontrol", "show", "job", j_id], capture_output=True, text=True, check=False)
            name = f"Job-{j_id}"
            for line in res.stdout.split():
                if line.startswith("JobName="):
                    name = line.split("=")[1]
                    break
            tracked_jobs[j_id] = name
            logger.info(f"Job submitted: {name} (ID: {j_id})")
        except Exception:
            tracked_jobs[j_id] = f"Job-{j_id}"
            logger.info(f"Job submitted: ID {j_id}")

    # Active monitoring tracking loop
    while tracked_jobs:
        time.sleep(30)
        try:
            job_str = ",".join(tracked_jobs.keys())
            # Fetch absolute status flags along with identifier groupings
            result = subprocess.run(
                ["squeue", "-j", job_str, "-h", "-o", "%A|%T"],
                capture_output=True, text=True, check=False
            )
            
            # Extract status mapping for items currently on the cluster matrix
            current_queue = {}
            for line in result.stdout.splitlines():
                if line.strip() and "|" in line:
                    j_id, state = line.strip().split("|", 1)
                    current_queue[j_id] = state

            # Find which jobs are no longer in the active queue map
            completed_this_round = []
            for active_id in list(tracked_jobs.keys()):
                if active_id not in current_queue:
                    # Double-check exit state via sacct to safely capture errors/crashes
                    exit_check = subprocess.run(
                        ["sacct", "-j", active_id, "-n", "-X", "--format=State"],
                        capture_output=True, text=True, check=False
                    )
                    final_state = exit_check.stdout.strip() if exit_check.stdout.strip() else "COMPLETED"
                    
                    if any(err in final_state for err in ["FAIL", "TIMEOUT", "NODE_FAIL", "CANCELLED"]):
                        logger.error(f"Pipeline stopped or error detected: {tracked_jobs[active_id]} (ID: {active_id}) ended with status {final_state}")
                        sys.exit(f"Pipeline execution halted due to cluster job failure on ID {active_id}.")
                    
                    logger.info(f"Job completed: {tracked_jobs[active_id]} (ID: {active_id})")
                    completed_this_round.append(active_id)

            # Drop completed IDs out of active checking matrix map
            for c_id in completed_this_round:
                del tracked_jobs[c_id]

        except Exception as e:
            logger.error(f"Error querying cluster status via squeue: {e}")
            time.sleep(10)


def read_sample_info(file_name):
    """Extracts Sample ID keys safely from the first column of the info file."""
    samples = []
    with open(file_name, 'r') as f:
        next(f)  # Skip header line
        for line in f:
            if line.strip():
                samples.append(line.split('\t')[0])
    return samples


def read_config(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)


def init_logging(log_path, file_name):
    logger = logging.getLogger('pipeline')
    logger.setLevel(logging.INFO)
    
    os.makedirs(log_path, exist_ok=True)
    fh = logging.FileHandler(os.path.join(log_path, file_name))
    sh = logging.StreamHandler()
    
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    sh.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(sh)

if __name__ == "__main__":
    main()