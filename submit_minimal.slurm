#!/bin/bash
#SBATCH --job-name=tcr_epitope_minimal_gpu
#SBATCH --output=/cluster/home/martihe1/log/%j_%N__tcr_minimal.out
#SBATCH --error=/cluster/home/martihe1/log/%j_%N__tcr_minimal.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --nodelist=sacramento
#SBATCH --account=cai_exp

echo "=========================================="
echo "Starting TCR-Epitope Minimal Baseline (GPU)"
echo "=========================================="
echo "Job ID:  $SLURM_JOB_ID"
echo "Node:    $(hostname)"
echo "Time:    $(date)"
echo "=========================================="

# ── environment ────────────────────────────────────────────────────────────
module load python/3.10.14
source /cluster/home/martihe1/venv/bin/activate

# install scirpy/shap if missing (lightweight check)
python3 -c "import scirpy" 2>/dev/null || pip install -q scirpy scanpy
python3 -c "import shap"   2>/dev/null || pip install -q shap

# ── output and working directory ───────────────────────────────────────────
cd /cluster/home/martihe1
mkdir -p /cluster/home/martihe1/outputs
mkdir -p /cluster/home/martihe1/log

# ── run ────────────────────────────────────────────────────────────────────
python3 -u /cluster/home/martihe1/codes/tcr_epitope_minimal.py

echo ""
echo "=========================================="
echo "Job finished at $(date)"
echo "=========================================="

# ── reminder to copy outputs to local ─────────────────────────────────────
echo ""
echo "To copy outputs to local machine run:"
echo "  scp -r martihe1@dallas.zhaw.ch:/cluster/home/martihe1/outputs/ ~/projectsZHAW/VT1/iggytop/outputs/"
