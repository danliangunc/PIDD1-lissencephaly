#!/bin/bash
#SBATCH --job-name=SubmitCellRangerJobs
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
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

filelist='3751D70Org1_CZ15_HHT 3751D70Org2_CZ16_HHT 3751D70Org3_CZ17_HHT Y6D70Org1_CZ1_HHT Y6D70Org2_CZ2_HHT Y6D70Org3_CZ3_HHT';

for file in $filelist
do
sbatch -o $file.out Run_CellRangerCount.sh $file;

done;

exit 0;
