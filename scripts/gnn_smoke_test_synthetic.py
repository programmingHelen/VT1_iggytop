"""
Smoke test for GNN training pipeline using synthetic data.
This bypasses the need for pre-computed ESM2 embeddings.

Usage:
  source .venv/bin/activate
  python3 scripts/gnn_smoke_test_synthetic.py
"""

import torch
import torch.nn as nn
import numpy as np
from torch_geometric.data import Data, DataLoader
from iggytop.gnn.models import build_model, GNNLinkPredictor
from sklearn.metrics import roc_auc_score

def create_synthetic_graph(num_tcr_nodes=100, num_epitope_nodes=50, num_edges=300, embedding_dim=320):
    """Create a synthetic bipartite TCR-Epitope graph."""
    num_nodes = num_tcr_nodes + num_epitope_nodes
    
    # Random node features
    x = torch.randn(num_nodes, embedding_dim, dtype=torch.float32)
    
    # Random bipartite edges (TCR -> Epitope)
    edge_list = []
    labels = []
    for _ in range(num_edges):
        tcr_idx = np.random.randint(0, num_tcr_nodes)
        epi_idx = num_tcr_nodes + np.random.randint(0, num_epitope_nodes)
        label = np.random.randint(0, 2)
        edge_list.append([tcr_idx, epi_idx])
        labels.append(label)
    
    edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
    y = torch.tensor(labels, dtype=torch.long)
    
    graph = Data(x=x, edge_index=edge_index, y=y, num_nodes=num_nodes)
    return graph


def test_training_loop():
    """Test a full training loop with synthetic data."""
    print("Creating synthetic graph...")
    graph = create_synthetic_graph(num_tcr_nodes=100, num_epitope_nodes=50, num_edges=300, embedding_dim=64)
    
    # Create loaders
    train_loader = DataLoader([graph], batch_size=32, shuffle=True)
    val_loader = DataLoader([graph], batch_size=32, shuffle=False)
    
    print(f"Graph: {graph.num_nodes} nodes, {graph.edge_index.shape[1]} edges")
    print()
    
    # Build model
    device = torch.device("cpu")
    model = build_model(
        "graphsage",
        input_dim=64,
        hidden_dim=32,
        output_dim=16,
        num_layers=2,
    ).to(device)
    
    predictor = GNNLinkPredictor(embedding_dim=16, hidden_dim=16).to(device)
    criterion = nn.BCELoss()
    optimizer = torch.optim.Adam(
        list(model.parameters()) + list(predictor.parameters()),
        lr=0.001
    )
    
    print(f"Model parameters: {sum(p.numel() for p in model.parameters() if p.requires_grad)}")
    print(f"Predictor parameters: {sum(p.numel() for p in predictor.parameters() if p.requires_grad)}")
    print()
    
    # Train one epoch
    print("Training for 1 epoch...")
    model.train()
    predictor.train()
    
    for batch_idx, batch in enumerate(train_loader):
        batch = batch.to(device)
        optimizer.zero_grad()
        
        # Forward
        node_emb = model(batch.x, batch.edge_index)
        
        # Simple link prediction: average TCR and epitope embeddings
        edge_list = batch.edge_index
        tcr_emb = node_emb[edge_list[0]]
        epi_emb = node_emb[edge_list[1]]
        
        logits = predictor(tcr_emb, epi_emb)
        labels = batch.y.unsqueeze(1).float()
        
        loss = criterion(logits, labels)
        loss.backward()
        optimizer.step()
        
        print(f"  Batch {batch_idx}: loss={loss.item():.4f}")
    
    # Validation
    print("\nValidating...")
    model.eval()
    predictor.eval()
    
    all_preds = []
    all_labels = []
    
    with torch.no_grad():
        for batch in val_loader:
            batch = batch.to(device)
            node_emb = model(batch.x, batch.edge_index)
            edge_list = batch.edge_index
            tcr_emb = node_emb[edge_list[0]]
            epi_emb = node_emb[edge_list[1]]
            logits = predictor(tcr_emb, epi_emb)
            all_preds.extend(logits.cpu().numpy().flatten())
            all_labels.extend(batch.y.numpy())
    
    all_preds = np.array(all_preds)
    all_labels = np.array(all_labels)
    
    auc = roc_auc_score(all_labels, all_preds)
    print(f"  Validation AUC: {auc:.4f}")
    print()
    print("✅ Training pipeline works!")


if __name__ == "__main__":
    test_training_loop()
