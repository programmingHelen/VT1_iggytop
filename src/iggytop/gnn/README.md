# GNN Module for TCR-Epitope Binding Prediction

A modular, easy-to-configure Graph Neural Network framework for link prediction in TCR-Epitope bipartite networks.

## Architecture Overview

```
src/iggytop/gnn/
├── models.py              # GraphSAGE, GAT, GCN implementations
├── train.py              # Training loop with MLflow integration
├── data_loader.py        # Graph construction and data loading
├── __init__.py
└── config/
    ├── graphsage.yaml    # GraphSAGE hyperparameters
    ├── gat.yaml          # GAT hyperparameters
    ├── gcn.yaml          # GCN hyperparameters
    └── train.yaml        # Training & data config
```

## Features

✅ **Three GNN architectures** with consistent interface
✅ **YAML configuration** for easy hyperparameter tuning
✅ **Early stopping** with configurable patience
✅ **MLflow integration** for experiment tracking
✅ **Bipartite graph** construction from AnnData
✅ **Batch normalization** and **residual connections** (GAT)
✅ **Link prediction head** for binding classification

## Installation

Add to `pyproject.toml`:

```toml
torch-geometric>=2.3.0
torch>=2.0.0
pyyaml>=6.0
mlflow>=2.0.0
```

Install:
```bash
pip install torch-geometric
pip install pyyaml mlflow
```

## Quick Start

### 1. Training with Command Line

```bash
# Train GraphSAGE with default config
python -m src.iggytop.gnn.train --config graphsage

# Train GAT with custom hyperparameters
python -m src.iggytop.gnn.train --config gat --epochs 100 --lr 0.001 --batch-size 64

# Train GCN on CPU
python -m src.iggytop.gnn.train --config gcn --device cpu

# Disable MLflow tracking
python -m src.iggytop.gnn.train --config graphsage --no-mlflow
```

### 2. Training Programmatically

```python
from iggytop.gnn import build_model, GNNTrainer, create_data_loaders
import yaml

# Load configs
with open("src/iggytop/gnn/config/graphsage.yaml") as f:
    model_config = yaml.safe_load(f)
with open("src/iggytop/gnn/config/train.yaml") as f:
    train_config = yaml.safe_load(f)

# Create data loaders
train_loader, val_loader, test_loader = create_data_loaders(
    anndata_path="data/deduplicated_anndata.h5ad",
    batch_size=32,
    seed=42
)

# Initialize trainer
trainer = GNNTrainer(model_config["model"], train_config)

# Train
model, predictor, history = trainer.train(
    train_loader, 
    val_loader, 
    test_loader,
    num_epochs=100
)
```

### 3. Building Models Directly

```python
from iggytop.gnn import build_model
import torch

# Build any model
model = build_model(
    "graphsage",
    input_dim=320,        # ESM2 embedding dimension
    hidden_dim=128,       # Hidden layer dimension
    output_dim=64,        # Output embedding dimension
    num_layers=3,         # Number of GNN layers
    aggregator_type="mean",
    dropout=0.3
)

# Forward pass
x = torch.randn(1000, 320)  # 1000 nodes with 320-dim features
edge_index = torch.randint(0, 1000, (2, 5000))  # 5000 edges
out = model(x, edge_index)  # [1000, 64] embeddings
```

## Configuration Files

### Model Config (graphsage.yaml, gat.yaml, gcn.yaml)

Key parameters:
- `model.num_layers`: Number of GNN layers (usually 2-4)
- `model.hidden_dim`: Hidden dimension (64-256)
- `model.output_dim`: Embedding dimension (32-128)
- `model.dropout`: Dropout rate (0.0-0.5)
- `training.learning_rate`: Adam learning rate
- `training.early_stopping_patience`: Epochs without improvement before stopping

Example modification (graphsage.yaml):
```yaml
model:
  num_layers: 4           # Add more layers for complex patterns
  hidden_dim: 256         # Increase capacity
  dropout: 0.2            # Reduce dropout for smaller datasets
training:
  learning_rate: 0.0005   # Lower learning rate
  early_stopping_patience: 20
```

### Train Config (train.yaml)

- `data.anndata_path`: Path to input AnnData file
- `data.embedding_type`: "esm2" for ESM2 embeddings
- `output.save_dir`: Where to save checkpoints
- `output.use_mlflow`: Enable/disable MLflow tracking

## Model Comparison

| Model      | Strengths                              | Best For                    |
|-----------|----------------------------------------|-----------------------------|
| **GraphSAGE** | Fast, scalable, inductive learning | Large graphs, new nodes     |
| **GAT**     | Interpretable attention, multi-head   | Understanding node importance |
| **GCN**     | Lightweight, classic spectral theory  | Smaller graphs, lower memory |

## Output

Checkpoints saved to `outputs/gnn_models/`:
- `{model}_best.pt` — Best model weights
- `{model}_history.json` — Training history (loss, AUC, etc.)

MLflow logged to:
- `outputs/mlflow.db` (configurable in train.yaml)
- View with: `mlflow ui --backend-store-uri sqlite:///outputs/mlflow.db`

## Hyperparameter Tuning

Start simple, then expand:

```yaml
# Conservative (likely to work)
model:
  num_layers: 2
  hidden_dim: 64
  dropout: 0.3
training:
  learning_rate: 0.001
  early_stopping_patience: 10

# Aggressive (for complex patterns)
model:
  num_layers: 4
  hidden_dim: 256
  dropout: 0.2
training:
  learning_rate: 0.0005
  early_stopping_patience: 20
```

Then run grid search:
```bash
for lr in 0.001 0.0005 0.0001; do
  for layers in 2 3 4; do
    python -m src.iggytop.gnn.train --config graphsage --lr $lr --epochs 100
  done
done
```

## API Reference

### `build_model(model_name, input_dim, hidden_dim, output_dim, num_layers, **kwargs)`
Factory function to instantiate any model.

### `BaseGNNModel`
Abstract base class. Subclassed by GraphSAGE, GAT, GCN.

- `forward(x, edge_index, batch=None)` → node embeddings

### `GNNTrainer`
Main trainer class.

- `train(train_loader, val_loader, test_loader, num_epochs)` → (model, predictor, history)
- `train_epoch(...)` → train metrics
- `evaluate(...)` → validation metrics

### `BipartiteGraphDataset`
Constructs graphs from AnnData files.

- `build_graph()` → PyTorch Geometric Data object
- `split_by_edges(graph, train_ratio, val_ratio)` → train/val/test graphs

### `GNNLinkPredictor`
Link prediction head.

- `forward(z1, z2)` → binding probabilities

## Citation

If you use this module, please cite:

- **GraphSAGE**: Hamilton et al., "Inductive Representation Learning on Large Graphs", NeurIPS 2017
- **GAT**: Veličković et al., "Graph Attention Networks", ICLR 2018
- **GCN**: Kipf & Welling, "Semi-Supervised Classification with Graph Convolutional Networks", ICLR 2017

## Troubleshooting

### Out of Memory
- Reduce `batch_size` in config
- Reduce `hidden_dim` or `num_layers`
- Use `device: cpu` for gradient accumulation

### Poor Convergence
- Lower `learning_rate`
- Increase `early_stopping_patience`
- Check data balance (positive/negative ratio)

### Slow Training
- Use GraphSAGE (faster than GAT/GCN)
- Reduce `num_layers`
- Enable multi-GPU with DataParallel (TODO)
