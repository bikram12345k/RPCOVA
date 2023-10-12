#!/bin/bash
#SBATCH --job-name=comparison_methods_ate_changesettings_revision_v4
#SBATCH --array=1-640
#SBATCH --partition=hpg2-compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL              # Mail events (NONE, BEGIN, END, FAIL, AL$
#SBATCH --mem=4gb				#Memory limit
#SBATCH --mail-user=bkarmakar@ufl.edu    # Where to send mail
#SBATCH --time=1:41:00               # Time limit hrs:min:sec
#SBATCH --account=bkarmakar
#SBATCH --qos=bkarmakar-b
##SBATCH --output=arrayid

date; hostname; pwd

export R_LIBS=/blue/bkarmakar/bkarmakar/struncatedP

module load R

Rscript --vanilla comparison_methods_ate_changesettings_revision_v4.R $SLURM_ARRAY_TASK_ID
