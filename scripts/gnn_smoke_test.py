"""
Lightweight smoke test to validate GNN model imports and forward passes.
Runs a synthetic forward pass for GraphSAGE, GAT, and GCN.

Usage:
  source .venv/bin/activate
  python3 scripts/gnn_smoke_test.py
"""

import torch
from iggytop.gnn.models import build_model

DEVICE = torch.device("cpu")

def test_model(model_name: str):
    print(f"Testing {model_name}...")
    model = build_model(
        model_name,
        input_dim=16,
        hidden_dim=8,
        output_dim=4,
        num_layers=2,
    )
    model.to(DEVICE)
    model.eval()

    # small synthetic graph: 10 nodes, 30 random edges
    x = torch.randn(10, 16, device=DEVICE)
    edge_index = torch.randint(0, 10, (2, 30), device=DEVICE)

    with torch.no_grad():
        out = model(x, edge_index)
    print(f"  output shape: {out.shape}\n")

if __name__ == "__main__":
    for name in ["graphsage", "gat", "gcn"]:
        try:
            test_model(name)
        except Exception as e:
            print(f"Error testing {name}: {e}")
            raise
    print("Smoke test finished.")
