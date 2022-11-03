#!/bin/bash

#SBATCH --output=../outputs/groups_%a.out
#SBATCH --nodes=1                     
#SBATCH --time=6:00:00   
#SBATCH --partition=high_priority
#SBATCH --account=gbesla  
#SBATCH --array 0-99
#SBATCH --qos=user_qos_gbesla

module load python/3.8

date
python3 utils/make_groups.py $SLURM_ARRAY_TASK_ID "TNG"
date

