#!/bin/bash

#SBATCH --output=../outputs/pair_med_%a.out
#SBATCH --nodes=1                     
#SBATCH --time=6:00:00   
#SBATCH --partition=high_priority
#SBATCH --account=gbesla  
#SBATCH --array 0-135
#SBATCH --qos=user_qos_gbesla

module load python/3.8

date
python3 utils/make_1pair.py $SLURM_ARRAY_TASK_ID "Illustris"

if [ $SLURM_ARRAY_TASK_ID -lt 100 ]
then
    echo 'thats less'
    python3 utils/make_1pair.py $SLURM_ARRAY_TASK_ID "TNG" 
fi 
date

