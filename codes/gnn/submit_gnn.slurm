#!/bin/bash
#SBATCH --job-name=tcr_gnn
#SBATCH --output=/cluster/home/martihe1/log/%j_%N__tcr_gnn.out
#SBATCH --error=/cluster/home/martihe1/log/%j_%N__tcr_gnn.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu_top
#SBATCH --account=cai_exp

echo "Job ID : $SLURM_JOB_ID"
echo "Node   : $(hostname)"
echo "Time   : $(date)"

module load python/3.10.14
source /cluster/home/martihe1/venv/bin/activate

# install torch_geometric if missing
python3 -c "import torch_geometric" 2>/dev/null || {
    echo "Installing torch_geometric ..."
    pip install -q torch_geometric
    pip install -q torch_scatter torch_sparse -f https://data.pyg.org/whl/torch-$(python3 -c "import torch; print(torch.__version__)").html
}

python3 -c "import yaml" 2>/dev/null || pip install -q pyyaml

cd /cluster/home/martihe1
mkdir -p outputs/checkpoints log

CONFIG=${1:-configs/sage_dot.yaml}   # default: sage_dot, override with: sbatch submit_gnn.slurm configs/gat_mlp.yaml
MODE=${2:-train}

echo "Config : $CONFIG"
echo "Mode   : $MODE"

# ── smoke test first ─────────────────────────────────────────────────────────
echo "Running smoke test ..."
python3 -u /cluster/home/martihe1/codes/gnn/train.py \
    --config /cluster/home/martihe1/codes/gnn/$CONFIG \
    --mode one_batch

if [ $? -ne 0 ]; then
    echo "Smoke test FAILED — aborting."
    exit 1
fi

echo "Smoke test passed. Starting training ..."

# ── full training ─────────────────────────────────────────────────────────────
python3 -u /cluster/home/martihe1/codes/gnn/train.py \
    --config /cluster/home/martihe1/codes/gnn/$CONFIG \
    --mode $MODE

echo "Done at $(date)"

# Submit examples:
#   sbatch submit_gnn.slurm configs/sage_dot.yaml train   ← fast baseline
#   sbatch submit_gnn.slurm configs/gat_mlp.yaml  train   ← attention model
#   sbatch submit_gnn.slurm configs/sage_dot.yaml training.epochs=5   ← quick test
