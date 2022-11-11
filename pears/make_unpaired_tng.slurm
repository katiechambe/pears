#!/bin/bash

#SBATCH --output=../outputs/pairs/unTNG_%a.out
#SBATCH --nodes=1                     
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00   
#SBATCH --partition=high_priority
#SBATCH --account=gbesla  
#SBATCH --qos=user_qos_gbesla
#SBATCH --array 0-96

module load python/3.8

date
python3 utils/make_unpaired.py $SLURM_ARRAY_TASK_ID "TNG" 10
date


