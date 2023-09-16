#!/bin/bash
#SBATCH --job-name=Submit_Seurat_Jobs
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

SampleIDs='3751D70Org1_CZ15 3751D70Org2_CZ16 Y5LOF12522-CZ23 C9Org12522-CZ21 3751D70Org3_CZ17 YB712522-CZ24 Y6D70Org1_CZ1 YB7HOMO12522-CZ22 YB7HomoOrg2-CZ33 Y6D70Org3_CZ3 Y6D70Org2_CZ2 C9HomoOrg1-CZ32';

#sbatch -o 3751D70Org1_CZ15.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R 3751D70Org1_CZ15";

sbatch -o 3751D70Org2_CZ16.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R 3751D70Org2_CZ16";

sbatch -o Y5LOF12522-CZ23.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R Y5LOF12522-CZ23";

sbatch -o C9Org12522-CZ21.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R C9Org12522-CZ21";

sbatch -o 3751D70Org3_CZ17.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R 3751D70Org3_CZ17";

sbatch -o YB712522-CZ24.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R YB712522-CZ24";

sbatch -o Y6D70Org1_CZ1.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R Y6D70Org1_CZ1";

sbatch -o YB7HOMO12522-CZ22.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R YB7HOMO12522-CZ22";

sbatch -o YB7HomoOrg2-CZ33.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R YB7HomoOrg2-CZ33";

sbatch -o Y6D70Org3_CZ3.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R Y6D70Org3_CZ3";

sbatch -o Y6D70Org2_CZ2.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R Y6D70Org2_CZ2";

sbatch -o C9HomoOrg1-CZ32.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R C9HomoOrg1-CZ32";

sbatch -o C9HomoOrg2.out -n 1 --mem=64g --time=2-00:00:00 --wrap="Rscript Seurat_snRNA.R C9HomoOrg2";


