import os
import yaml
import subprocess

def run_step(workdir, sample_info, cfg_file, samples):
    cfg = yaml.safe_load(open(cfg_file))
    BAM_DIR = os.path.join(workdir, "bam")
    TMP = cfg["paths"]["temp"]
    FG, PICARD, BWA, BWA_INDEX, SAMTOOLS, REF = cfg["tools"]["fgbio"], cfg["tools"]["picard"], cfg["tools"]["bwa"], cfg["references"]["bwa_index"], cfg["tools"]["samtools"], cfg["references"]["fasta"]

    submitted_jobs = []
    for s in samples:
        job = f"{workdir}/job_files/group.by.umi.{s}.sh"
        with open(job, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f"#SBATCH --job-name=umi.{s}\n")
            f.write(f"#SBATCH --output={workdir}/job_output/group.by.umi.{s}.%j.out\n")
            f.write(f"#SBATCH --error={workdir}/job_error/group.by.umi.{s}.%j.err\n")
            f.write("#SBATCH --ntasks=12\n\nset -euo pipefail\n\n")
            
            t = f"{TMP}/group_reads_by_umi/{s}"
            f.write(f"rm -rf {t} && mkdir -p {t}/bam\n")
            f.write(f"cp {BAM_DIR}/{s}.mapped.zipperbam.bam {t}/bam/input.bam\nln -sf {t}/bam/input.bam {t}/bam/{s}.sort.bam\n")
            f.write(f"java -Xmx64g -jar {PICARD} CollectHsMetrics I={t}/bam/{s}.sort.bam O={t}/bam/{s}_hs_before.txt R={REF} BAIT_INTERVALS={cfg['panels']['onco_interval']} TARGET_INTERVALS={cfg['panels']['onco_interval']}\n")
            f.write(f"cp {t}/bam/{s}_hs_before.txt {BAM_DIR}/\n")
            f.write(f"java -jar {FG} GroupReadsByUmi --input {t}/bam/{s}.sort.bam --strategy Adjacency --edits 5 --output {t}/bam/{s}.grouped.bam\n")
            f.write(f"java -jar {FG} CallMolecularConsensusReads --input {t}/bam/{s}.grouped.bam --output {t}/bam/{s}.consensus.unmapped.bam --min-reads 1 --min-input-base-quality 20 --threads 24\n")
            f.write(f"samtools fastq {t}/bam/{s}.consensus.unmapped.bam | {BWA} mem -t 24 -K 150000000 -p -Y {BWA_INDEX} - | samtools view -b -o {t}/bam/{s}.consensus.bam -\n")
            f.write(f"samtools sort -n --no-PG -@ 24 -o {t}/bam/{s}.consensus.qname.bam {t}/bam/{s}.consensus.bam\n")
            f.write(f"samtools sort -n --no-PG -@ 24 -o {t}/bam/{s}.unmapped.qname.bam {t}/bam/{s}.consensus.unmapped.bam\n")
            f.write(f"java -jar {FG} ZipperBams --input {t}/bam/{s}.consensus.qname.bam --unmapped {t}/bam/{s}.unmapped.qname.bam --ref {REF} --output {t}/bam/{s}.consensus.mapped.bam\n")
            f.write(f"java -jar {FG} FilterConsensusReads --input {t}/bam/{s}.consensus.mapped.bam --output /dev/stdout --ref {REF} --min-reads 1 --min-base-quality 20 --max-base-error-rate 0.2 | samtools sort --threads 24 -o {t}/bam/{s}.consensus.mapped.filtered.bam\n")
            f.write(f"cp {t}/bam/{s}.consensus.mapped.filtered.bam {BAM_DIR}/\n")
            f.write(f"java -Xmx64g -jar {PICARD} CollectHsMetrics -I {t}/bam/{s}.consensus.mapped.filtered.bam -O {t}/bam/{s}_hs_after.txt -R {REF} --BAIT_INTERVALS {cfg['panels']['onco_interval']} --TARGET_INTERVALS {cfg['panels']['onco_interval']} --PER_TARGET_COVERAGE {t}/bam/{s}_per_target_coverage.txt --COMPRESSION_LEVEL 1 --MAX_RECORDS_IN_RAM 50000000\n")
            f.write(f"cp {t}/bam/{s}_hs_after.txt {BAM_DIR}/\n")
            f.write(f"cp {t}/bam/{s}_per_target_coverage.txt {BAM_DIR}/\n")
            f.write(f"rm -rf {t}\n")

        os.chmod(job, 0o755)
        res = subprocess.run(["sbatch", "--parsable", job], capture_output=True, text=True)
        if res.returncode == 0:
            submitted_jobs.append(res.stdout.strip())
            
    return submitted_jobs