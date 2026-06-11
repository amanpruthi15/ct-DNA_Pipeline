import os
import glob
import yaml
import subprocess


def run_step(workdir, sample_info, config_file, samples):

    with open(config_file) as f:
        cfg = yaml.safe_load(f)

    TRIMMOMATIC = cfg["tools"]["trimmomatic"]

    ADAPTER = cfg.get(
        "references", {}
    ).get(
        "adapters",
        "/home/act/database/adapters.fa"
    )

    trim = cfg.get(
        "trimming",
        {
            "illumina_clip": {
                "seed_mismatches": 2,
                "palindrome_clip_threshold": 30,
                "simple_clip_threshold": 10,
                "min_adapter_length": 8,
                "keep_both_reads": True
            },
            "headcrop": 0,
            "trailing": 3,
            "slidingwindow": {
                "window": 4,
                "quality": 15
            },
            "minlen": 36
        }
    )

    ILLUMINACLIP = (
        f"ILLUMINACLIP:{ADAPTER}:"
        f"{trim['illumina_clip']['seed_mismatches']}:"
        f"{trim['illumina_clip']['palindrome_clip_threshold']}:"
        f"{trim['illumina_clip']['simple_clip_threshold']}:"
        f"{trim['illumina_clip']['min_adapter_length']}:"
        f"{str(trim['illumina_clip']['keep_both_reads']).lower()}"
    )

    HEADCROP = trim["headcrop"]
    TRAILING = trim["trailing"]
    SW = trim["slidingwindow"]["window"]
    Q = trim["slidingwindow"]["quality"]
    MINLEN = trim["minlen"]

    FASTQ_DIR = os.path.join(workdir, "fastq")
    TRIM_DIR = os.path.join(workdir, "trimmed_fastq")

    JOB_DIR = os.path.join(workdir, "job_files")
    OUT_DIR = os.path.join(workdir, "job_output")
    ERR_DIR = os.path.join(workdir, "job_error")

    os.makedirs(TRIM_DIR, exist_ok=True)
    os.makedirs(JOB_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(ERR_DIR, exist_ok=True)

    submitted_jobs = []

    def write_trim_command(fq1, fq2, prefix):
        return (
            f"{TRIMMOMATIC} PE "
            f"{fq1} {fq2} "
            f"{prefix}.R1.paired.fastq.gz "
            f"{prefix}.R1.unpaired.fastq.gz "
            f"{prefix}.R2.paired.fastq.gz "
            f"{prefix}.R2.unpaired.fastq.gz "
            f"{ILLUMINACLIP} "
            f"HEADCROP:{HEADCROP} "
            f"TRAILING:{TRAILING} "
            f"SLIDINGWINDOW:{SW}:{Q} "
            f"MINLEN:{MINLEN}\n"
        )

    for sample in samples:

        r1_files = sorted(
            glob.glob(
                os.path.join(
                    FASTQ_DIR,
                    f"{sample}*R1*.fastq.gz"
                )
            )
        )

        r2_files = sorted(
            glob.glob(
                os.path.join(
                    FASTQ_DIR,
                    f"{sample}*R2*.fastq.gz"
                )
            )
        )

        if len(r1_files) == 0:
            print(f"[SKIP] {sample}: no FASTQ files found")
            continue

        if len(r1_files) != len(r2_files):
            print(
                f"[SKIP] {sample}: "
                f"R1 count ({len(r1_files)}) != "
                f"R2 count ({len(r2_files)})"
            )
            continue

        print(
            f"[INFO] {sample}: "
            f"found {len(r1_files)} lane(s)"
        )

        job_file = os.path.join(
            JOB_DIR,
            f"trim.{sample}.sh"
        )

        with open(job_file, "w") as out:

            out.write("#!/bin/bash\n")
            out.write(f"#SBATCH --job-name=trim.{sample}\n")
            out.write(f"#SBATCH --output={OUT_DIR}/trim.{sample}.%j.out\n")
            out.write(f"#SBATCH --error={ERR_DIR}/trim.{sample}.%j.err\n")
            out.write(f"#SBATCH --ntasks=4\n")
            out.write("\n")

            out.write("set -euo pipefail\n\n")

            out.write(f"cd {TRIM_DIR}\n\n")

            if len(r1_files) == 1:

                out.write(
                    write_trim_command(
                        r1_files[0],
                        r2_files[0],
                        sample
                    )
                )

                out.write(
                    f"mv {sample}.R1.paired.fastq.gz "
                    f"{sample}.trimmed.merged_R1.fastq.gz\n"
                )

                out.write(
                    f"mv {sample}.R2.paired.fastq.gz "
                    f"{sample}.trimmed.merged_R2.fastq.gz\n"
                )

            else:

                lane_prefixes = []

                for idx, (fq1, fq2) in enumerate(
                    zip(r1_files, r2_files),
                    start=1
                ):

                    lane_prefix = f"{sample}.L{idx}"
                    lane_prefixes.append(lane_prefix)

                    out.write(
                        write_trim_command(
                            fq1,
                            fq2,
                            lane_prefix
                        )
                    )

                merged_r1 = " ".join(
                    [
                        f"{p}.R1.paired.fastq.gz"
                        for p in lane_prefixes
                    ]
                )

                merged_r2 = " ".join(
                    [
                        f"{p}.R2.paired.fastq.gz"
                        for p in lane_prefixes
                    ]
                )

                out.write(
                    f"\ncat {merged_r1} > "
                    f"{sample}.trimmed.merged_R1.fastq.gz\n"
                )

                out.write(
                    f"cat {merged_r2} > "
                    f"{sample}.trimmed.merged_R2.fastq.gz\n"
                )

                out.write(
                    f"rm -f {' '.join([f'{p}.R1.paired.fastq.gz' for p in lane_prefixes])}\n"
                )

                out.write(
                    f"rm -f {' '.join([f'{p}.R2.paired.fastq.gz' for p in lane_prefixes])}\n"
                )

            out.write(
                f"rm -f {sample}*.unpaired*.fastq.gz\n"
            )

        os.chmod(job_file, 0o755)

        res = subprocess.run(
            ["sbatch", "--parsable", job_file],
            capture_output=True,
            text=True
        )

        if res.returncode == 0:

            job_id = res.stdout.strip()

            print(
                f"[SUBMITTED] {sample}: "
                f"job {job_id}"
            )

            submitted_jobs.append(job_id)

        else:

            print(
                f"[FAILED] {sample}: "
                f"{res.stderr.strip()}"
            )

    return submitted_jobs