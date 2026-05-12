"""
Graph Neural Network module for TCR-Epitope binding prediction.

Implements:
- GraphSAGE: Inductive representation learning
- GAT: Graph Attention Networks
- GCN: Graph Convolutional Networks

Usage:
    from iggytop.gnn import build_model, GNNTrainer
    
    # Build a model
    model = build_model("graphsage", input_dim=320, hidden_dim=128, output_dim=64, num_layers=3)
    
    # Or use the trainer
    trainer = GNNTrainer(model_config, train_config)
    model, predictor, history = trainer.train(train_loader, val_loader, test_loader)
"""

from .models import GraphSAGE, GAT, GCN, GNNLinkPredictor, build_model, BaseGNNModel
from .train import GNNTrainer
from .data_loader import BipartiteGraphDataset, create_data_loaders

__all__ = [
    "GraphSAGE",
    "GAT",
    "GCN",
    "GNNLinkPredictor",
    "build_model",
    "BaseGNNModel",
    "GNNTrainer",
    "BipartiteGraphDataset",
    "create_data_loaders",
]

__version__ = "0.1.0"
