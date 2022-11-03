#!/bin/bash

#SBATCH --output=../outputs/ammass_ill_%a.out
#SBATCH --nodes=1                     
#SBATCH --time=6:00:00   
#SBATCH --partition=high_priority
#SBATCH --account=gbesla  
#SBATCH --array 0-135
#SBATCH --qos=user_qos_gbesla

module load python/3.8

date
python3 utils/make_am_realizations.py $SLURM_ARRAY_TASK_ID "Illustris" "dark" 1000
python3 utils/make_am_realizations.py $SLURM_ARRAY_TASK_ID "Illustris" "hydro" 1000
date

