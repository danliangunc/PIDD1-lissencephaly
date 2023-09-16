#!/bin/bash
#SBATCH --job-name=Run_CellRangerCount
#SBATCH --time=5-00:00:00
#SBATCH --ntasks=2 --cpus-per-task=16
#SBATCH --mem=112G
#SBATCH --mail-type=FAIL

mem_bytes=$(</sys/fs/cgroup/memory/slurm/uid_${SLURM_JOB_UID}/job_${SLURM_JOB_ID}/memory.limit_in_bytes)
mem_gbytes=$(( $mem_bytes / 1024 **3 ))

echo "Starting at $(date)"
echo "Job submitted to the ${SLURM_JOB_PARTITION} partition"
echo "Job name: ${SLURM_JOB_NAME}, Job ID: ${SLURM_JOB_ID}"
echo "  I have ${SLURM_CPUS_ON_NODE} CPUs and ${mem_gbytes}GiB of RAM on compute node $(hostname)"

## Run cellranger count to align the reads to hg38 ref genome
## not use "--include-introns" to calculate reads mapped to introns
## to keep cellranger as the same setting with sequencing core

ID=$1

/home/dl2248/project/Apps/cellranger-6.1.2/cellranger count \
	--id=$ID\_cellranger \
	--fastqs=/home/dl2248/project/CoCeZhang/FASTQ \
	--sample=$ID \
	--transcriptome=/home/dl2248/project/CoCeZhang/CellRanger/refdata-gex-GRCh38-2020-A

exit 0;

