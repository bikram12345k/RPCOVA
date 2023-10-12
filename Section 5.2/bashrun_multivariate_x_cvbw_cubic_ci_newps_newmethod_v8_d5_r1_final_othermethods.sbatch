#!/bin/bash
#SBATCH --job-name=RPCOVA_ci_simulation_cte_standard_error_ci_cubic_newps_newmethod_v8_d5_r1_fixll_fixedbw_othermethods
#SBATCH --array=1-390
#SBATCH --partition=hpg2-compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL              # Mail events (NONE, BEGIN, END, FAIL, AL$
#SBATCH --mem=7gb				#Memory limit
#SBATCH --mail-user=bkarmakar@ufl.edu    # Where to send mail
#SBATCH --time=30:41:00               # Time limit hrs:min:sec
#SBATCH --account=bkarmakar
#SBATCH --qos=bkarmakar-b
##SBATCH --output=arrayid

date; hostname; pwd

export PATH="${PATH}:/home/bkarmakar/Rpersonal/bin"
export R_LIBS=/blue/bkarmakar/bkarmakar/struncatedP
module load gurobi

FILE=simulation_cte_standard_error_ci_cubic.R
if test -f "$FILE"; then
    echo "$FILE exists."
fi

R --slave --no-restore --file=simulation_cte_standard_error_ci_cubic_newps_newmethod_v8_d5_r1_fixll_draft_final.R --args $SLURM_ARRAY_TASK_ID
