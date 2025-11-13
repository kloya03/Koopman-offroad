#!/bin/bash

#SBATCH --job-name=Sl_1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=56
#SBATCH --mem=350G
#SBATCH --time=30:00:00

module add matlab/2023b

cd /scratch/kloya/Koopman-offroad/scripts/simulation/palmetto_files

matlab -nodisplay -nosplash < batch1_sandyloam.m > results_sandyloam_1.txt