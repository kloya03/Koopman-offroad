#!/bin/bash
#SBATCH --job-name=test_array
#SBATCH --array=1-3                     # 3 tasks (same as number of lines in params_test.txt)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=test_output_%A_%a.txt  # Unique output per task

module add matlab/2025a

cd $SLURM_SUBMIT_DIR

PARAM_FILE="params_test.txt"

# Extract the line matching this array index
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PARAM_FILE")

# Parse a and b
read a b <<< "$LINE"

echo "Task ${SLURM_ARRAY_TASK_ID} reading line: $LINE"

# Create tag
TAG="a${a}_b${b}"

# Run MATLAB with the parameters
matlab -nodisplay -nosplash -r "a=${a}; b=${b}; param_tag=sprintf('%s','${TAG}'); test_array; exit;"

