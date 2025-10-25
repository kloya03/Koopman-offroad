#!/bin/bash

#SBATCH --job-name=Sl_3
#SBATCH --nodes=1
#SBATCH --tasks-per-node=56
#SBATCH --mem=350G
#SBATCH --time=30:00:00

module add matlab/2023b

cd /scratch/kloya/Koopman-offroad/scripts/simulation

matlab -nodisplay -nosplash < batch3_sandyloam.m > results_sandyloam_3.txt