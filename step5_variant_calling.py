import os
import yaml
import subprocess

def run_step(project_dir, sample_info_file, config_file, samples):
    with open(config_file) as f:
        cfg = yaml.safe_load(f)

    # Establish structural subdirectories
    job_files_dir = os.path.join(project_dir, "job_files")
    job_output_dir = os.path.join(project_dir, "job_output")
    job_error_dir = os.path.join(project_dir, "job_error")
    mutect2_dir = os.path.join(project_dir, "mutect2")
    freebayes_dir = os.path.join(project_dir, "freebayes")
    merged_variants_dir = os.path.join(project_dir, "merged_variants")

    for d in [mutect2_dir, freebayes_dir, merged_variants_dir]:
        os.makedirs(d, exist_ok=True)

    gatk_path = cfg["tools"]["gatk"]
    freebayes_path = cfg["tools"]["freebayes"]
    bcftools_path = cfg["tools"]["bcftools"]
    bgzip_path = cfg["tools"]["bgzip"]
    tabix_path = cfg["tools"]["tabix"]
    reference_fasta = cfg["references"]["fasta"]
    onco_bed = cfg["panels"]["onco_bed"]
    onco_interval = cfg["panels"]["onco_interval"]
    pon = cfg["mutect2"]["panel_of_normals"]
    germline = cfg["mutect2"]["germline_resource"]
    funcotator_ds = cfg["funcotator"]["datasource"]
    vep_cache = cfg["vep"]["cache_dir"]
    vep_plugins = cfg["vep"]["plugin_dir"]

    mutect2_job_ids = []
    freebayes_job_ids = []
    final_funcotator_job_ids = []

    # Submit Mutect2
    for sample in samples:
        job_file_path = os.path.join(job_files_dir, f"{sample}_mutect2.sh")
        with open(job_file_path, 'w') as f:
            f.write(f"#!/bin/bash\n#SBATCH --job-name=mutect2.{sample}\n#SBATCH --output={job_output_dir}/mutect2.{sample}.%j.out\n#SBATCH --error={job_error_dir}/mutect2.{sample}.%j.err\n#SBATCH --ntasks=2\n\n")
            f.write(f"samtools addreplacerg -r ID:{sample} -r SM:{sample} -o {project_dir}/bam/{sample}.consensus.recal.filt.rg.bam {project_dir}/bam/{sample}.consensus.recal.filt.bam\n")
            f.write(f"samtools index {project_dir}/bam/{sample}.consensus.recal.filt.rg.bam\n")
            f.write(f"{gatk_path} Mutect2 -R {reference_fasta} -I {project_dir}/bam/{sample}.consensus.recal.filt.rg.bam -L {onco_interval} --tumor-lod-to-emit 0.0 --minimum-allele-fraction 0.0 --f1r2-tar-gz {mutect2_dir}/{sample}.f1r2.tar.gz --panel-of-normals {pon} --germline-resource {germline} -O {mutect2_dir}/{sample}.mutect2.vcf.gz\n")
            f.write(f"{gatk_path} LearnReadOrientationModel -I {mutect2_dir}/{sample}.f1r2.tar.gz -O {mutect2_dir}/{sample}.artifact-priors.tar.gz\n")
            f.write(f"{gatk_path} FilterMutectCalls -R {reference_fasta} -V {mutect2_dir}/{sample}.mutect2.vcf.gz --orientation-bias-artifact-priors {mutect2_dir}/{sample}.artifact-priors.tar.gz -O {mutect2_dir}/{sample}.mutect2.filtered.vcf.gz\n")
            f.write(f"rm -f {project_dir}/bam/{sample}.consensus.recal.filt.rg.bam*\n")

        os.chmod(job_file_path, 0o755)
        res = subprocess.run(f"sbatch --parsable {job_file_path}", shell=True, capture_output=True, text=True)
        if res.returncode == 0:
            mutect2_job_ids.append(res.stdout.strip())

    # Build Dependency String for FreeBayes/VC Concat
    m2_deps = ":".join(mutect2_job_ids)

    for sample in samples:
        fb_job_path = os.path.join(job_files_dir, f"{sample}_freebayes.sh")
        with open(fb_job_path, 'w') as f:
            f.write(f"#!/bin/bash\n#SBATCH --job-name=freebayes.{sample}\n#SBATCH --output={job_output_dir}/freeabayes.{sample}.%j.out\n#SBATCH --error={job_error_dir}/freebayes.{sample}.%j.err\n#SBATCH --ntasks=2\n\n")
            f.write(f"{freebayes_path} --fasta-reference {reference_fasta} --bam {project_dir}/bam/{sample}.consensus.recal.filt.bam --targets {onco_bed} --min-alternate-fraction 0.0001 --min-alternate-count 5 --pooled-discrete --ploidy 2 > {freebayes_dir}/{sample}.freebayes.vcf\n")
            f.write(f"{bcftools_path} sort -o {freebayes_dir}/{sample}.freebayes.sorted.vcf {freebayes_dir}/{sample}.freebayes.vcf\n")
            f.write(f"{bgzip_path} {freebayes_dir}/{sample}.freebayes.sorted.vcf\n")
            f.write(f"{tabix_path} {freebayes_dir}/{sample}.freebayes.sorted.vcf.gz\n")
            f.write(f"{bcftools_path} concat -a {mutect2_dir}/{sample}.mutect2.filtered.vcf.gz {freebayes_dir}/{sample}.freebayes.sorted.vcf.gz -Oz -o {merged_variants_dir}/{sample}.merged.vcf.gz\n")
            f.write(f"{bcftools_path} norm -d none {merged_variants_dir}/{sample}.merged.vcf.gz -o {merged_variants_dir}/{sample}.normalized.vcf\n")

        os.chmod(fb_job_path, 0o755)
        res = subprocess.run(f"sbatch --dependency=afterok:{m2_deps} --parsable {fb_job_path}", shell=True, capture_output=True, text=True)
        if res.returncode == 0:
            freebayes_job_ids.append(res.stdout.strip())

    fb_deps = ":".join(freebayes_job_ids)

    # Submit Downstream Annotations (VEP followed by Funcotator)
    for sample in samples:
        vep_job_path = os.path.join(job_files_dir, f"{sample}_vep.sh")
        with open(vep_job_path, 'w') as f:
            f.write(f"#!/bin/bash\n#SBATCH --job-name=vep.{sample}\n#SBATCH --output={job_output_dir}/vep.{sample}.%j.out\n#SBATCH --error={job_error_dir}/vep.{sample}.%j.err\n#SBATCH --ntasks=2\n\nsource ~/.bashrc\nconda activate vep_env\n")
            f.write(f"vep --input_file {merged_variants_dir}/{sample}.normalized.vcf --output_file {merged_variants_dir}/variants_{sample}.vep.vcf --format vcf --vcf --symbol --terms SO --tsl --biotype --hgvs --species homo_sapiens --fasta {reference_fasta} --offline --cache --dir_cache {vep_cache} --dir_plugins {vep_plugins} --plugin ReferenceQuality --everything --force_overwrite\n")

        os.chmod(vep_job_path, 0o755)
        vep_res = subprocess.run(f"sbatch --dependency=afterok:{fb_deps} --parsable {vep_job_path}", shell=True, capture_output=True, text=True)
        
        if vep_res.returncode == 0:
            vep_id = vep_res.stdout.strip()
            
            func_job_path = os.path.join(job_files_dir, f"{sample}_funcotator.sh")
            with open(func_job_path, 'w') as f:
                f.write(f"#!/bin/bash\n#SBATCH --job-name=func_{sample}\n#SBATCH --output={job_output_dir}/funcotator.{sample}.%j.out\n#SBATCH --error={job_error_dir}/funcotator.{sample}.%j.err\n#SBATCH --ntasks=2\n\n")
                f.write(f"{gatk_path} Funcotator --variant {merged_variants_dir}/{sample}.normalized.vcf --reference {reference_fasta} --ref-version hg38 --data-sources-path {funcotator_ds} --output {merged_variants_dir}/{sample}.funcotated.vcf --output-file-format VCF --disable-sequence-dictionary-validation\n")

            os.chmod(func_job_path, 0o755)
            func_res = subprocess.run(f"sbatch --dependency=afterok:{vep_id} --parsable {func_job_path}", shell=True, capture_output=True, text=True)
            if func_res.returncode == 0:
                final_funcotator_job_ids.append(func_res.stdout.strip())

    # Return the terminal layer IDs to the orchestrator checking loop
    return final_funcotator_job_ids