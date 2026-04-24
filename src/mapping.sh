#!/bin/bash
#SBATCH --job-name=rnaseq_nfcore
#SBATCH --output=/home/allenma/Mary_Allen/voles2025/eando/slurm-%A.out
#SBATCH --error=/home/allenma/Mary_Allen/voles2025/eando/slurm-%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --cpus-per-task=16                  # Number of cores for main Nextflow job (not individual alignment)
#SBATCH --mem=40G                           # 32–40GB for main Nextflow manager job
#SBATCH --time=24:00:00                     # 24h is typical for pipelines
#SBATCH --qos=normal

module load nextflow/24.04.4

resultsdir=/pl/active/DonaldsonLab/Mary_Allen/voles2025/nfcoresrnaseq10092025/
samplesheet=${resultsdir}samplesheet.csv
gtf=/pl/active/DonaldsonLab/Mary_Allen/voles2025/reference/Microtus_ochrogaster.MicOch1.0.115.gtf
fasta=/pl/active/DonaldsonLab/Mary_Allen/voles2025/reference/Microtus_ochrogaster.MicOch1.0.dna.toplevel.fa
cf=nextflow.config


mkdir -p $SLURM_SCRATCH/nextflow/

nextflow run nf-core/rnaseq -r 3.17.0 --input $samplesheet --outdir $SLURM_SCRATCH/nextflow/ --fasta $fasta --gtf $gtf -profile apptainer -c $cf

rsync -r $SLURM_SCRATCH/nextflow/ $resultsdir
