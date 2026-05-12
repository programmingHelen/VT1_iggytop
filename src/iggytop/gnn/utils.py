"""
Utility functions for GNN training and evaluation.
"""

import torch
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Optional
import logging

logger = logging.getLogger(__name__)


def set_seed(seed: int):
    """Set random seed for reproducibility."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    logger.info(f"Random seed set to {seed}")


def count_parameters(model: torch.nn.Module) -> int:
    """Count total number of parameters in model."""
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def save_checkpoint(
    epoch: int,
    model: torch.nn.Module,
    predictor: torch.nn.Module,
    optimizer: torch.optim.Optimizer,
    metrics: Dict,
    save_path: str,
):
    """Save training checkpoint."""
    checkpoint = {
        "epoch": epoch,
        "model_state": model.state_dict(),
        "predictor_state": predictor.state_dict(),
        "optimizer_state": optimizer.state_dict(),
        "metrics": metrics,
    }
    torch.save(checkpoint, save_path)
    logger.info(f"Checkpoint saved to {save_path}")


def load_checkpoint(
    checkpoint_path: str,
    model: torch.nn.Module,
    predictor: torch.nn.Module,
    optimizer: Optional[torch.optim.Optimizer] = None,
    device: str = "cpu",
) -> Dict:
    """Load training checkpoint."""
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint["model_state"])
    predictor.load_state_dict(checkpoint["predictor_state"])
    if optimizer:
        optimizer.load_state_dict(checkpoint["optimizer_state"])
    logger.info(f"Checkpoint loaded from {checkpoint_path}")
    return checkpoint


def compute_metrics_summary(history: Dict) -> pd.DataFrame:
    """Convert training history to summary DataFrame."""
    rows = []
    for epoch, (train_metrics, val_metrics) in enumerate(zip(history["train"], history["val"])):
        row = {"epoch": epoch + 1}
        row.update({f"train_{k}": v for k, v in train_metrics.items()})
        row.update({f"val_{k}": v for k, v in val_metrics.items()})
        rows.append(row)
    return pd.DataFrame(rows)


def get_best_epoch(history: Dict, metric: str = "val_auc") -> int:
    """Find epoch with best validation metric."""
    val_metrics = history["val"]
    scores = [m.get(metric.replace("val_", ""), 0) for m in val_metrics]
    best_idx = np.argmax(scores)
    return best_idx + 1, scores[best_idx]


def create_node_embeddings_df(
    model: torch.nn.Module,
    graph_data: "Data",
    node_types: np.ndarray,
    device: str = "cpu",
) -> pd.DataFrame:
    """
    Generate embeddings for all nodes and return as DataFrame.
    
    Args:
        model: Trained GNN model
        graph_data: PyTorch Geometric graph
        node_types: Array indicating node type (0=TCR, 1=Epitope)
        device: Device to use
    
    Returns:
        DataFrame with columns: node_id, node_type, embedding
    """
    model.eval()
    with torch.no_grad():
        x = graph_data.x.to(device)
        edge_index = graph_data.edge_index.to(device)
        embeddings = model(x, edge_index).cpu().numpy()
    
    node_type_names = np.array(["TCR", "Epitope"])[node_types]
    
    df = pd.DataFrame({
        "node_id": np.arange(len(embeddings)),
        "node_type": node_type_names,
        **{f"emb_{i}": embeddings[:, i] for i in range(embeddings.shape[1])}
    })
    
    return df


def get_device() -> torch.device:
    """Get appropriate device (GPU if available, else CPU)."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"Using device: {device}")
    return device


def format_metrics(metrics: Dict, decimal_places: int = 4) -> str:
    """Format metrics dictionary as readable string."""
    lines = []
    for key, value in sorted(metrics.items()):
        if isinstance(value, float):
            lines.append(f"  {key:<20}: {value:.{decimal_places}f}")
        else:
            lines.append(f"  {key:<20}: {value}")
    return "\n".join(lines)


def get_top_important_nodes(
    model: torch.nn.Module,
    graph_data: "Data",
    node_type: str = "TCR",
    k: int = 10,
    device: str = "cpu",
) -> pd.DataFrame:
    """
    Get top-K nodes by embedding magnitude (rough importance proxy).
    
    Note: This is a heuristic. True importance requires model-specific
    interpretation (e.g., SHAP values, attention weights for GAT).
    """
    model.eval()
    with torch.no_grad():
        x = graph_data.x.to(device)
        edge_index = graph_data.edge_index.to(device)
        embeddings = model(x, edge_index).cpu().numpy()
    
    magnitudes = np.linalg.norm(embeddings, axis=1)
    
    node_type_mask = graph_data.node_type == (0 if node_type == "TCR" else 1)
    node_type_indices = np.where(node_type_mask.numpy())[0]
    
    top_idx = np.argsort(magnitudes[node_type_indices])[-k:][::-1]
    actual_idx = node_type_indices[top_idx]
    
    return pd.DataFrame({
        "node_id": actual_idx,
        "magnitude": magnitudes[actual_idx],
        "embedding_dim": [embeddings.shape[1]] * len(actual_idx),
    })
