#!/bin/bash
#SBATCH --job-name=gnn_one_epoch
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --account=YOUR_ACCOUNT_HERE

set -e

# Edit these paths as needed for your cluster user
PROJECT_DIR="$HOME/projectsZHAW/VT1/iggytop"
VENV_PATH="$HOME/.venv"  # or point to your cluster venv
OUTPUT_DIR="$PROJECT_DIR/outputs/gnn_slurm"

mkdir -p "$OUTPUT_DIR"
cd "$PROJECT_DIR"

# Activate venv
if [ -f "$VENV_PATH/bin/activate" ]; then
    source "$VENV_PATH/bin/activate"
else
    echo "Virtualenv not found at $VENV_PATH. Activate environment manually or edit the script." >&2
    exit 1
fi

echo "Running on $(hostname)"
echo "PWD=$(pwd)"

# Optional: install requirements in cluster venv if missing (commented out by default)
# pip install -r requirements_cluster.txt

python3 -m src.iggytop.gnn.train --config graphsage --epochs 1 --no-mlflow

# Copy or leave outputs in $OUTPUT_DIR as needed
