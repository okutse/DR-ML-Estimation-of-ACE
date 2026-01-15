#!/bin/bash
#SBATCH --job-name=autogluon_stage
#SBATCH --output=logs/autogluon_stage_%j.out
#SBATCH --error=logs/autogluon_stage_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=$USER@jhu.edu

# Exit on error
set -e

# Print job info
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "========================================"

# Change to project root directory
cd /users/iyliu/Code/DR-ML-Estimation-of-ACE

# Create logs directory if it doesn't exist
mkdir -p logs

# Activate environment (adjust as needed for your setup)
# Option 1: Using uv
source .venv/bin/activate

# Option 2: Using conda (uncomment if using conda)
# source ~/.bashrc
# conda activate autogluon_env

which python

# Run the training script
python scripts/train_autogluon_stage.py \
    --data data/merged_data.csv \
    --sample data/data_clinical_sample.txt \
    --output results/autogluon_stage \
    --num-folds 10 \
    --presets good_quality \
    --seed 42 \
    --secondary-stratify CANCER_TYPE_DETAILED ER_STATUS \
    --threshold-metric f1 \
    --max-shap-samples 200 \
    --max-lime-samples 5 \
    --permutation-subsample 2000

echo "========================================"
echo "End time: $(date)"
echo "========================================"
