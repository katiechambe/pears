#!/bin/bash

#SBATCH --output=../outputs/pairs/tng1K_%a.out
#SBATCH --nodes=1                     
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=2
#SBATCH --time=60:00:00   
#SBATCH --partition=high_priority
#SBATCH --account=gbesla  
#SBATCH --qos=user_qos_gbesla
#SBATCH --array 0-99

module load python/3.8

date
echo "starting with pairs"
python3 utils/make_pairs.py $SLURM_ARRAY_TASK_ID "TNG" 1000
echo "finished collecting pairs"
echo "starting with unpaired"
python3 utils/make_unpaired.py $SLURM_ARRAY_TASK_ID "TNG" 1000
echo "finished with unpaired"
date


