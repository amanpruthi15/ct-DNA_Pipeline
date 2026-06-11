import os
import yaml
import subprocess

def run_step(workdir, sample_info, cfg_file, samples):
    cfg = yaml.safe_load(open(cfg_file))
    BAM_DIR = os.path.join(workdir, "bam")
    GATK, REF, DBSNP = cfg["tools"]["gatk"], cfg["references"]["fasta"], cfg["known_sites"]["dbsnp"]

    submitted_jobs = []
    for s in samples:
        job = f"{workdir}/job_files/base.recaliberation.{s}.sh"
        with open(job, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f"#SBATCH --job-name=bqsr.{s}\n")
            f.write(f"#SBATCH --output={workdir}/job_output/base.recaliberation.{s}.%j.out\n")
            f.write(f"#SBATCH --error={workdir}/job_error/base.recaliberation.{s}.%j.err\n")
            f.write("#SBATCH --ntasks=8\n\nset -euo pipefail\n\n")
            
            f.write(f"samtools sort --no-PG -@ 24 -o {BAM_DIR}/{s}.consensus.mapped.filtered.coord.bam {BAM_DIR}/{s}.consensus.mapped.filtered.bam\n")
            f.write(f"samtools index {BAM_DIR}/{s}.consensus.mapped.filtered.coord.bam\n")
            f.write(f"{GATK} BaseRecalibrator -I {BAM_DIR}/{s}.consensus.mapped.filtered.coord.bam -R {REF} --known-sites {DBSNP} -O {BAM_DIR}/{s}.recal.table\n")
            f.write(f"{GATK} ApplyBQSR -R {REF} -I {BAM_DIR}/{s}.consensus.mapped.filtered.coord.bam --bqsr-recal-file {BAM_DIR}/{s}.recal.table -O {BAM_DIR}/{s}.recal.bam\n")
            f.write(f"samtools view -b -L {cfg['panels']['onco_bed']} -o {BAM_DIR}/{s}.consensus.recal.filt.bam {BAM_DIR}/{s}.recal.bam\n")
            f.write(f"samtools index {BAM_DIR}/{s}.consensus.recal.filt.bam\n")
            f.write(f"samtools index {BAM_DIR}/{s}.consensus.recal.filt.bam\n")
            f.write(f"rm -f {BAM_DIR}/{s}.recal.bam {BAM_DIR}/{s}.recal.bai {BAM_DIR}/{s}.consensus.mapped.filtered.coord.bam*\n")

        os.chmod(job, 0o755)
        res = subprocess.run(["sbatch", "--parsable", job], capture_output=True, text=True)
        if res.returncode == 0:
            submitted_jobs.append(res.stdout.strip())
            
    return submitted_jobs