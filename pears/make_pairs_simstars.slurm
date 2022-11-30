#!/bin/bash

#SBATCH --output=../outputs/pairs/simstars_%a.out
#SBATCH --nodes=1                     
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00   
#SBATCH --partition=high_priority
#SBATCH --account=gbesla  
#SBATCH --qos=user_qos_gbesla
#SBATCH --array 0-135

module load python/3.8

date
echo "starting Illustris"
###python3 utils/make_pairs_simstars.py $SLURM_ARRAY_TASK_ID "Illustris"
echo "finished collecting halos in Illustris"


if [ $SLURM_ARRAY_TASK_ID -lt 100 ]
then
    echo "starting TNG"
    python3 utils/make_pairs_simstars.py $SLURM_ARRAY_TASK_ID "TNG"
    echo "finished with TNG"
fi 

date


