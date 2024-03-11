#!/bin/bash

#Written by: libr8020@colorado.edu
# Date: 20231114
# Purpose: for doing shuffled gene correlations

#SBATCH --partition=amilan     # Alpine partition
#SBATCH --qos=normal                 # Summit qos
#SBATCH --time=00:03:00           # Max wall time
#SBATCH --ntasks=1           # Number of tasks per job
#SBATCH --nodes=1          # nodes
#SBATCH --mem=1G
#SBATCH --job-name=run_correlations        # Job submission name
#SBATCH --output=/projects/libr8020/run_R_scripts/eo_files/run_correlations.%j.out   # Output file name with Job ID
#SBATCH --error=/projects/libr8020/run_R_scripts/eo_files/run_correlations.%j.err
#SBATCH --mail-type=ALL            # Email user when job finishes
#SBATCH --mail-user=libr8020@colorado.edu # Email address of user

# purge all existing modules
module purge

# load conda environment
source /curc/sw/anaconda/default
conda activate r_env

#cd /projects/libr8020/run_R_scripts/outdir/
#outdir=/projects/libr8020/run_R_scripts/outdir/
#start R
Rscript --vanilla correlation_oncluster_withingrps.R $clust $gene $outdir


