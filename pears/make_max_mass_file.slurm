#!/bin/bash

#SBATCH --output=../outputs/tng-maxmass%a.out
#SBATCH --nodes=1                     
#SBATCH --time=10:00:00   
#SBATCH --partition=high_priority
#SBATCH --account=gbesla  
#SBATCH --array 0-135
#SBATCH --qos=user_qos_gbesla

module load python/3.8

date
python3 utils/make_max_mass_file.py $SLURM_ARRAY_TASK_ID Illustris
date

