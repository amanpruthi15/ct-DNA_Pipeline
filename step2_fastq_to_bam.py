import os
import glob
import yaml
import subprocess
from collections import defaultdict

def run_step(workdir, sample_info, config_file):
    with open(config_file) as f:
        cfg = yaml.safe_load(f)

    FASTQ_DIR = os.path.join(workdir, cfg["paths"]["fastq"])
    BAM_DIR = os.path.join(workdir, cfg["paths"]["bam"])
    TMP_ROOT = cfg["paths"]["temp"]
    FG, BWA, SAMTOOLS, REF, BWA_INDEX = cfg["tools"]["fgbio"], cfg["tools"]["bwa"], cfg["tools"]["samtools"], cfg["references"]["fasta"], cfg["references"]["bwa_index"]

    fq_dict = defaultdict(dict)
    for fq in glob.glob(f"{FASTQ_DIR}/*.fastq.gz"):
        base = os.path.basename(fq)
        if "_R1.fastq.gz" in base:
            sample = base.split(".trimmed")[0]
            fq_dict[sample]["R1"] = fq
        elif "_R2.fastq.gz" in base:
            sample = base.split(".trimmed")[0]
            fq_dict[sample]["R2"] = fq

    submitted_jobs = []
    for sample, reads in fq_dict.items():
        if "R1" not in reads or "R2" not in reads:
            continue

        r1, r2 = reads["R1"], reads["R2"]
        job = f"{workdir}/job_files/fastq2bam.{sample}.sh"

        with open(job, "w") as out:
            out.write("#!/bin/bash\n")
            out.write(f"#SBATCH --job-name=fastq2bam.{sample}\n")
            out.write(f"#SBATCH --output={workdir}/job_output/fastq2bam.{sample}.%j.out\n")
            out.write(f"#SBATCH --error={workdir}/job_error/fastq2bam.{sample}.%j.err\n")
            out.write("#SBATCH --nodes=1\n#SBATCH --ntasks=12\n\n")
            out.write("set -euo pipefail\n\n")
            
            tmp = f"{TMP_ROOT}/fastq_to_bam/{sample}"
            out.write(f'TMP={tmp}\nmkdir -p $TMP/bam\nmkdir -p {BAM_DIR}\ntrap "rm -rf $TMP" EXIT\n')
            out.write(f'java -jar {FG} FastqToBam --input {r1} {r2} --read-structures 5M2S+T 5M2S+T --sample {sample} --library {sample} --output $TMP/bam/{sample}.unmapped.bam\n')
            out.write(f'cp $TMP/bam/{sample}.unmapped.bam {BAM_DIR}/\n')
            out.write(f'samtools fastq $TMP/bam/{sample}.unmapped.bam | {BWA} mem -t 24 -K 150000000 -p -Y {BWA_INDEX} - | samtools view -b -o $TMP/bam/{sample}.mapped.bam -\n')
            out.write(f'cp $TMP/bam/{sample}.mapped.bam {BAM_DIR}/\n')
            out.write(f'samtools sort --no-PG -n -@ 24 -o $TMP/bam/{sample}.mapped.qname.bam $TMP/bam/{sample}.mapped.bam\n')
            out.write(f'samtools sort --no-PG -n -@ 24 -o $TMP/bam/{sample}.unmapped.qname.bam $TMP/bam/{sample}.unmapped.bam\n')
            out.write(f'java -jar {FG} ZipperBams --input $TMP/bam/{sample}.mapped.qname.bam --unmapped $TMP/bam/{sample}.unmapped.qname.bam --ref {REF} --output $TMP/bam/{sample}.mapped.zipperbam.bam\n')
            out.write(f'cp $TMP/bam/{sample}.mapped.zipperbam.bam {BAM_DIR}/\n')

        os.chmod(job, 0o755)
        res = subprocess.run(["sbatch", "--parsable", job], capture_output=True, text=True)
        if res.returncode == 0:
            submitted_jobs.append(res.stdout.strip())
            
    return submitted_jobs