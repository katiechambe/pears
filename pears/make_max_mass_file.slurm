#!/bin/bash

#SBATCH --output=../outputs/maxmass_%a.out
#SBATCH --nodes=1                     
#SBATCH --time=6:00:00   
#SBATCH --partition=standard
#SBATCH --account=gbesla  
#SBATCH --array 130-135
#SBATCH --qos=user_qos_gbesla

module load python/3.8

date
python3 utils/make_max_mass_file.py $SLURM_ARRAY_TASK_ID Illustris > howmany.txt
date

