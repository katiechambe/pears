#!/bin/bash

#SBATCH --output=../outputs/makedata.out
#SBATCH --nodes=1                     
#SBATCH --time=10:00:00   
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --partition=high_priority
#SBATCH --account=gbesla 
#SBATCH --qos=user_qos_gbesla

module load python/3.8

date
python3 /xdisk/gbesla/katiechambe/pears/notebooks/makedata_slurm.py
 
date

