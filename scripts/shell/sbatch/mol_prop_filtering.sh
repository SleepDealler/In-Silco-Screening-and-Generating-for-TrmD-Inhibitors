#!/bin/bash

#SBATCH --partition=common

#SBATCH --qos=8gpu3d

#SBATCH --nodes=1

#SBATCH --ntasks=4

#SBATCH --gres=gpu:0

#SBATCH --job-name=ipz

#SBATCH --time=20:00:00

#SBATCH --output=/home/pszmkptrzk/run_logs/outputs/%x-%j.out

#SBATCH --error=/home/pszmkptrzk/run_logs/errors/%x-%j.error

#SBATCH --export=ALL

# Append CUDA path to the existing PATH
export PATH="/usr/local/cuda/bin:${PATH}"

# ----------------------------------------------------------------------------------------
SINGULARITY_CONTAINER="/home/pszmkptrzk/apptainer/initial-screening-container.sif"
# ----------------------------------------------------------------------------------------


# Set script name
SCRIPT_NAME=$(basename "$0")
CURRENT_DATE=$(date +%Y-%m-%d)

# Print communicates
echo "Script name ${SCRIPT_NAME}"
echo "Current date ${CURRENT_DATE}"
echo "Home directory: ${HOME}"
echo "Working directory: $PWD"
echo "Current node: ${SLURM_NODELIST}"
echo "PATH: $PATH"

# Print job information
echo "Current node: ${SLURM_NODELIST}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "GPUs allocated: ${SLURM_GPUS_PER_NODE}"

# Compute  info
echo "Number of CPU cores available: $(nproc)"
echo "Listing available GPUs with nvidia-smi:"
nvidia-smi
echo "Compact GPU list:"
nvidia-smi -L
echo "CUDA version:"
nvcc -V

# Run bash script in the singularity container
srun singularity exec --nv $SINGULARITY_CONTAINER bash -c "
    python3 --version && \
    conda list && \
    python3 src/run_molecular_property_based_filtering.py
"
# srun singularity exec --nv $SINGULARITY_CONTAINER bash -c "conda list"
# ". /opt/conda/etc/profile.d/conda.sh && conda activate ipz && conda list | grep rdkit"