# ctDNA Liquid Biopsy Variant Calling Pipeline

An automated, orchestrator-driven genomic analysis pipeline developed for high-throughput identification of low-frequency somatic variants from circulating tumor DNA (ctDNA) next-generation sequencing (NGS) data. The pipeline handles raw data trimming, unique molecular identifier (UMI) consensus tracking, alignment, base quality score recalibration, somatic variant calling, and functional annotation on a SLURM-managed High-Performance Computing (HPC) cluster.

## Pipeline Architecture & Features

The pipeline is split into 5 modular execution layers coordinated by a centralized runtime orchestrator (`pipeline_runner.py`):

1. **Step 1: Read Trimming (`step1_trimming.py`)** Uses Trimmomatic for adapter clipping, base-quality filtering, and sliding-window read pruning. Supports automatic lane merging for dual-lane datasets.
2. **Step 2: Fastq to BAM Alignment (`step2_fastq_to_bam.py`)** Converts raw FASTQ records into unmapped BAM files while preserving UMI structures via `fgbio`, followed by initial alignment using `bwa mem`.
3. **Step 3: UMI Molecular Consensus (`step3_umi_align_consensus.py`)** Groups reads by UMI adjacency constraints, calculates molecular duplex/single-strand consensus reads using `fgbio CallMolecularConsensusReads`, and remaps consensus sequences back to the reference genome to eliminate PCR and sequencing artifacts.
4. **Step 4: Base Quality Score Recalibration (`step4_base_recaliberation.py`)** Performs GATK BQSR on the coordinate-sorted consensus BAMs using known polymorphic sites (dbSNP) to minimize systemic base-calling bias.
5. **Step 5: Somatic Variant Calling & Annotation (`step5_variant_calling.py`)** Executes somatic single-nucleotide variant (SNV) and small indel discovery across Mutect2 and FreeBayes. Final variant files are cleanly normalized and functionally annotated through Ensembl VEP and GATK Funcotator.

### Features
* **Smart Live Monitoring:** Tracks job execution dynamically through `squeue` and `scontrol`. Logs are printed cleanly *only* when a job is submitted, completed, or errors out.
* **Granular Stage Resuming:** Each pipeline stage is gated behind a boolean flag in a centralized YAML file. You can enable, disable, or restart from any specific step effortlessly.
* **Fault-Tolerant Halting:** If a job fails or exits with a non-zero cluster state (`FAILED`, `TIMEOUT`, `NODE_FAIL`), the orchestrator captures it via `sacct` and halts down-stream execution instantly to protect computation resources.

---

## Directory Structure

Set up your project workspace matching the architecture expected by the wrapper scripts:

```text
/path/to/workdir/
├── fastq/                             # Place raw Illumina fastq.gz pairs here
│   ├── CF0201110_S4_R1_001.fastq.gz
│   └── CF0201110_S4_R2_001.fastq.gz
├── config.yaml                        # Pipeline parameter and tool path config
├── 2024-11-04_cfDNA_Val2_Sample.txt   # Tab-delimited sample tracker index file
├── ctDNA_pipeline_runner.py           # Core Master Orchestrator
├── step1_trimming.py
├── step2_fastq_to_bam.py
├── step3_umi_align_consensus.py
├── step4_base_recaliberation.py
└── step5_variant_calling.py
