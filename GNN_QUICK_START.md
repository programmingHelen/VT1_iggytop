# GNN Module — Quick Start Guide

## 📁 What Was Created

```
src/iggytop/gnn/
├── __init__.py              # Package initialization
├── models.py                # GraphSAGE, GAT, GCN implementations
├── train.py                 # Training loop with MLflow tracking
├── data_loader.py           # Graph construction from AnnData
├── utils.py                 # Helper functions
├── README.md                # Full documentation
└── config/
    ├── graphsage.yaml       # GraphSAGE hyperparameters
    ├── gat.yaml             # GAT hyperparameters
    ├── gcn.yaml             # GCN hyperparameters
    └── train.yaml           # Training & data configuration
```

## 🚀 Getting Started

### Step 1: Install Dependencies

Add to your `pyproject.toml`:

```toml
torch>=2.0.0
torch-geometric>=2.3.0
pyyaml>=6.0
mlflow>=2.0.0
scikit-learn>=1.0.0
```

Install:
```bash
pip install torch torch-geometric pyyaml mlflow
```

### Step 2: Run Training

**Option A: Command Line (Simplest)**
```bash
# Train GraphSAGE
cd ~/projectsZHAW/VT1/iggytop
python -m src.iggytop.gnn.train --config graphsage

# Train GAT with custom params
python -m src.iggytop.gnn.train --config gat --epochs 100 --lr 0.001

# Train GCN on CPU
python -m src.iggytop.gnn.train --config gcn --device cpu
```

**Option B: Python Script**
```python
from iggytop.gnn import GNNTrainer, create_data_loaders
import yaml

# Load configs
with open("src/iggytop/gnn/config/graphsage.yaml") as f:
    model_config = yaml.safe_load(f)
with open("src/iggytop/gnn/config/train.yaml") as f:
    train_config = yaml.safe_load(f)

# Create loaders
train_loader, val_loader, test_loader = create_data_loaders(
    anndata_path="data/deduplicated_anndata.h5ad",
    batch_size=32
)

# Train
trainer = GNNTrainer(model_config["model"], train_config)
model, predictor, history = trainer.train(train_loader, val_loader, test_loader)
```

### Step 3: Hyperparameter Tuning

Edit the YAML files to adjust parameters:

```yaml
# src/iggytop/gnn/config/graphsage.yaml
model:
  num_layers: 3              # ← Change this
  hidden_dim: 128            # ← Or this
  output_dim: 64
  dropout: 0.3               # ← Or this
  
training:
  learning_rate: 0.001       # ← Or this
  early_stopping_patience: 10
```

Then retrain to see if performance improves.

### Step 4: View Results

**Check MLflow Dashboard:**
```bash
mlflow ui --backend-store-uri sqlite:///outputs/mlflow.db
```

**View Saved Checkpoints:**
```bash
ls outputs/gnn_models/
# graphsage_best.pt, gat_best.pt, gcn_best.pt
# graphsage_history.json, etc.
```

## 🔧 Architecture Comparison

| Aspect | GraphSAGE | GAT | GCN |
|--------|-----------|-----|-----|
| Speed | ⚡ Fast | 🐢 Slower | ⚡ Fast |
| Memory | 💾 Low | 💾 Higher | 💾 Low |
| Interpretability | 🤔 Medium | 👁️ High (attention weights) | 🤔 Medium |
| Best for | Large graphs | Small precise graphs | Spectral understanding |

**Recommendation:** Start with **GraphSAGE** (fastest, then tune).

## 📊 Key Config Parameters Explained

### Model Configuration

```yaml
model:
  num_layers: 3              # Depth of network (2-4 usually sufficient)
  hidden_dim: 128            # Capacity (bigger = more expressive, slower)
  output_dim: 64             # Embedding dimension (node representation size)
  dropout: 0.3               # Prevent overfitting (0.2-0.4 good range)
  batch_norm: true           # Stabilize training
```

### Training Configuration

```yaml
training:
  batch_size: 32             # Smaller = more frequent updates, noisier gradient
  num_epochs: 100            # Max training steps
  learning_rate: 0.001       # Update step size (0.0001-0.01 typical)
  weight_decay: 1e-5         # L2 regularization
  early_stopping_patience: 10 # Stop if no improvement for N epochs
  val_split: 0.2             # Use 20% of data for validation
```

## 💡 Tips for Best Results

### For Small Datasets (< 1000 samples)
```yaml
model:
  num_layers: 2
  hidden_dim: 64
  dropout: 0.3
training:
  learning_rate: 0.001
  early_stopping_patience: 15
```

### For Large Datasets (> 10000 samples)
```yaml
model:
  num_layers: 4
  hidden_dim: 256
  dropout: 0.2
training:
  learning_rate: 0.0005
  batch_size: 128
  early_stopping_patience: 20
```

### If Overfitting
- Increase `dropout` (0.4-0.5)
- Decrease `hidden_dim`
- Add more `weight_decay`
- Use fewer `num_layers`

### If Underfitting
- Decrease `dropout`
- Increase `hidden_dim`
- Increase `num_layers`
- Decrease `learning_rate` (allows more careful training)

## 🧪 Experimentation Workflow

1. **Start simple:**
   ```bash
   python -m src.iggytop.gnn.train --config graphsage --epochs 50
   ```

2. **Check metrics** in `outputs/gnn_models/graphsage_history.json`

3. **If AUC < 0.55:** Increase model capacity
   - Edit `graphsage.yaml`: `hidden_dim: 256`

4. **If AUC plateaus early:** Longer training
   - Edit `graphsage.yaml`: `early_stopping_patience: 20`

5. **If training is slow:** Use GraphSAGE (if not already)
   - Switch to `--config graphsage`

6. **Once stable, fine-tune:**
   ```bash
   python -m src.iggytop.gnn.train --config gat --lr 0.0005 --epochs 150
   ```

## 📈 Expected Performance

From your analysis notebook, RF achieved **AUC ≈ 0.738** on concatenated ESM2 embeddings.

GNN expectations:
- **GraphSAGE**: AUC 0.70-0.75 (competitive with RF)
- **GAT**: AUC 0.72-0.77 (with proper tuning)
- **GCN**: AUC 0.65-0.72 (more conservative)

The advantage: GNNs can **leverage graph structure** (which RF cannot), potentially discovering new TCR-Epitope relationships.

## 🆘 Common Issues

**Issue: "CUDA out of memory"**
- Reduce `batch_size` (32 → 16)
- Reduce `hidden_dim` (128 → 64)
- Use `--device cpu`

**Issue: Training too slow**
- Use GraphSAGE instead of GAT
- Reduce `num_layers` (3 → 2)
- Increase `batch_size` (32 → 64)

**Issue: Poor convergence**
- Lower `learning_rate` (0.001 → 0.0005)
- Increase `early_stopping_patience` (10 → 20)
- Check data balance

**Issue: Model not loading**
- Make sure `torch-geometric` is installed correctly
- Try: `pip install --upgrade torch-geometric`

## 📚 Next Steps

1. ✅ Train baseline models (GraphSAGE, GAT, GCN)
2. 🎯 Compare with RF from `full_analysis.ipynb`
3. 🔍 Analyze attention weights (GAT) for interpretability
4. 🌳 Try GCN with pre-computed adjacency matrix
5. 🧩 Combine with other modalities (MHC class, species)

## Questions?

- Check `README.md` for detailed API documentation
- Look at `config/*.yaml` for all available parameters
- Run `python -m src.iggytop.gnn.train --help` for CLI options

Good luck with training! 🚀
