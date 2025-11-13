#!/bin/bash

#SBATCH --job-name=Cl_2
#SBATCH --nodes=1
#SBATCH --tasks-per-node=56
#SBATCH --mem=350G
#SBATCH --time=30:00:00

module add matlab/2023b

cd /scratch/kloya/Koopman-offroad/scripts/simulation/palmetto_files

matlab -nodisplay -nosplash < batch2_clay.m > results_clay_2.txt