#!/bin/bash
#SBATCH --job-name=cov_selection
#SBATCH --output=logs/covariate_selection_%j.out
#SBATCH --error=logs/covariate_selection_%j.err
#SBATCH --time=06:00:00
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
echo "Job Name: $SLURM_JOB_NAME (Full Covariate Selection)"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "========================================"

# Change to project root directory
cd /users/mliang23/DR-ML-Estimation-of-ACE

# Create logs directory if it doesn't exist
mkdir -p logs

# Activate virtual environment
echo "Activating virtual environment..."
source .venv/bin/activate

# Verify Python and packages
echo "Python location: $(which python)"
echo "Python version: $(python --version)"

# Check if required packages are installed
echo "Checking required packages..."
python -c "import autogluon.tabular; print('✓ AutoGluon installed')" || echo "WARNING: AutoGluon not found"
python -c "import sklearn; print('✓ sklearn installed')" || echo "WARNING: sklearn not found"
python -c "import shap; print('✓ SHAP installed')" || echo "WARNING: SHAP not found"
python -c "import lime; print('✓ LIME installed')" || echo "WARNING: LIME not found"

# Run full covariate selection with refit_full strategy
echo "========================================"
echo "Starting Full Covariate Selection..."
echo "  Strategy: AutoGluon with refit_full()"
echo "  - Train/Test split: 80/20"
echo "  - AutoGluon internal validation for model selection"
echo "  - refit_full() on 100% of training data"
echo "  - Feature importance from refitted model"
echo "  - SHAP analysis for covariate selection"
echo "  - Final evaluation on hold-out test set"
echo "  - Time limit: 3600s (1 hour)"
echo "  - Preset: good_quality"
echo "========================================"

python scripts/train_covariate_autogluon.py \
    --data-dir data \
    --output-dir results/covariate_selection_full \
    --use-nested-cv \
    --test-size 0.2 \
    --presets good_quality \
    --time-limit 3600 \
    --seed 42 \
    --max-shap-samples 200 \
    --max-lime-samples 5 \
    --eval-metric roc_auc

echo "========================================"
echo "Job completed: $(date)"
echo "========================================"

# Print results
if [ -f results/covariate_selection_full/final_results.json ]; then
    echo ""
    echo "📊 Final Results:"
    cat results/covariate_selection_full/final_results.json
fi

# Print selected features
if [ -f results/covariate_selection_full/selected_features.json ]; then
    echo ""
    echo "✅ Selected Features:"
    cat results/covariate_selection_full/selected_features.json
fi

echo ""
echo "📁 All results saved to: results/covariate_selection_full/"
