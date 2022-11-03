#!/bin/bash

#SBATCH --output=../outputs/subs/ill_%a.out
###SBATCH --output=../outputs/subs/tng_%a.out
#SBATCH --nodes=1                     
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00   
#SBATCH --partition=high_priority
#SBATCH --account=gbesla  
#SBATCH --array 0-135
#SBATCH --qos=user_qos_gbesla

module load python/3.8

date
python3 utils/make_subhalos.py $SLURM_ARRAY_TASK_ID "Illustris" 
date

