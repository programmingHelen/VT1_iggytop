#!/usr/bin/env bash
set -euo pipefail

# Activate project venv and run a 1-epoch training to smoke-test the GNN pipeline.
# Usage: ./scripts/run_gnn_one_epoch.sh

# Adjust this if your venv lives elsewhere
VENV_PATH=".venv"
if [ ! -d "$VENV_PATH" ]; then
  echo "Virtualenv not found at $VENV_PATH. Activate your environment manually and run the command below instead:" >&2
  echo "  python3 -m src.iggytop.gnn.train --config graphsage --epochs 1 --no-mlflow" >&2
  exit 1
fi

# Source the venv
source "$VENV_PATH/bin/activate"

# Ensure we're running from project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR/.."
cd "$PROJECT_ROOT"

echo "Running one-epoch GNN training (GraphSAGE) from: $(pwd)"
python3 -m src.iggytop.gnn.train --config graphsage --epochs 1 --no-mlflow
