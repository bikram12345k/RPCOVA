#!/bin/bash
#SBATCH --job-name=multivariate_x_cvbw_linear_red_stress_test_d_revision_v4_rerun
#SBATCH --array=1-640
#SBATCH --partition=hpg2-compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL              # Mail events (NONE, BEGIN, END, FAIL, AL$
#SBATCH --mem=8gb				#Memory limit
#SBATCH --mail-user=bkarmakar@ufl.edu    # Where to send mail
#SBATCH --time=90:41:00               # Time limit hrs:min:sec
#SBATCH --account=bkarmakar
#SBATCH --qos=bkarmakar-b
##SBATCH --output=arrayid

date; hostname; pwd

export R_LIBS=/blue/bkarmakar/bkarmakar/struncatedP
export PATH="${PATH}:/home/bkarmakar/Rpersonal/bin"

module load gurobi


Rscript --vanilla simulation_multivariate_x_cvbw_linear_red_stress_test_d_revision_v4_bartps.R $SLURM_ARRAY_TASK_ID
