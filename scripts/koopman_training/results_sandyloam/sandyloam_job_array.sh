#!/bin/bash
#SBATCH --job-name=train_SL
#SBATCH --array=1-120                      # One task per parameter set
#SBATCH --nodes=1
#SBATCH --tasks-per-node=28
#SBATCH --mem=100G
#SBATCH --time=30:00:00
#SBATCH --output=part1/slurm_outputs/%x_%A_%a.txt            # Saves output in the folder you call sbatch from

module add matlab/2023b

# The directory where you call sbatch from:
SLURM_DIR="$SLURM_SUBMIT_DIR"


# MATLAB script directory:
MATLAB_DIR="$SLURM_DIR/.."
PARAM_FILE="$MATLAB_DIR/params.txt"
DATA_FILE="$MATLAB_DIR/../../datasets/sandyloam_100hz_no_elev_experiment_1579.mat"
FUNCTIONS_FILE="$MATLAB_DIR/../../functions/utility"

# Move to MATLAB script directory so it can find the .m file properly:
cd "$MATLAB_DIR"

# Extract the correct line (strip CR)
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PARAM_FILE" | tr -d '\r')

# Parse parameters: nB nl sy cut_off
read nl sy nB cut_off <<< "$LINE"

echo "Running Task ${SLURM_ARRAY_TASK_ID} with:"
echo "nB=$nB, nl=$nl, sy=$sy, cut_off=$cut_off"
echo "Using parameter file: $PARAM_FILE"

# Create tag
TAG="nB${nB}_nl${nl}_sy${sy}_c${cut_off}"

# Pass parameters AND the full param file path into MATLAB
matlab -nodisplay -nosplash -r "\
    nB=${nB}; \
    nl=${nl}; \
    sy=${sy}; \
    cut_off=${cut_off}; \
    param_tag=sprintf('%s','${TAG}'); \
    data_file_path=sprintf('%s','${DATA_FILE}'); \
    function_file_path=sprintf('%s','${FUNCTIONS_FILE}'); \
    Offroad_Koopman_Training_main; \
    exit;"
