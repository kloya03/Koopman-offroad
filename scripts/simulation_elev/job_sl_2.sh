#!/bin/bash

#SBATCH --job-name=Sl_2_elev
#SBATCH --nodes=1
#SBATCH --tasks-per-node=56
#SBATCH --mem=350G
#SBATCH --time=30:00:00

module add matlab/2023b

cd /scratch/kloya/Koopman-offroad/scripts/simulation_elev

matlab -nodisplay -nosplash < batch2_sandyloam.m > results_sandyloam_2.txt