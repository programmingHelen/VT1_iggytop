#!/bin/bash
#SBATCH --job-name=tcr_epitope_minimal
#SBATCH --output=/cluster/home/martihe1/outputs/slurm_%j.log
#SBATCH --error=/cluster/home/martihe1/outputs/slurm_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:a100sxm:1
#SBATCH --partition=gpu
#SBATCH --account=cai_exp

set -e

mkdir -p /cluster/home/martihe1/outputs

cd /cluster/home/martihe1

source /cluster/home/martihe1/venv/bin/activate

echo "Running on $(hostname)"
echo "PWD=$(pwd)"

if ! python3 -c "import scanpy" &>/dev/null; then
    echo "Installing dependencies..."
    pip install -r /cluster/home/martihe1/requirements_cluster.txt
fi

python3 codes/tcr_epitope_minimal.py

# ── reminder to copy outputs to local ─────────────────────────────────────
echo ""
echo "To copy outputs to local machine run:"
echo "  scp -r martihe1@dallas.zhaw.ch:/cluster/home/martihe1/outputs/ ~/projectsZHAW/VT1/iggytop/outputs/"
