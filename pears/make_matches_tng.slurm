#!/bin/bash

#SBATCH --output=../outputs/matches/TNG100_%a.out
#SBATCH --nodes=1                     
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00   
#SBATCH --partition=high_priority
#SBATCH --account=gbesla  
#SBATCH --qos=user_qos_gbesla
#SBATCH --array 0-99

module load python/3.8

date
python3 utils/make_matched_pairs.py $SLURM_ARRAY_TASK_ID "TNG" 100
date


