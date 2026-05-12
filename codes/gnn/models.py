"""
models.py
=========
GNN models for TCR–Epitope binding prediction (link prediction on bipartite graph).

Architecture
------------
  Bipartite graph:
    - TCR nodes   (features = ESM2 CDR3β embedding, 320-dim)
    - Epitope nodes (features = ESM2 epitope embedding, 320-dim)
    - Edges = known binding pairs (positive examples from training set)

  Two encoder options:
    - GraphSAGE: inductive, works on unseen nodes at test time
    - GAT: attention-weighted aggregation over neighbours

  Link prediction head:
    - Dot product (fast, parameter-free)
    - MLP (learnable, more expressive)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv, GATConv, to_hetero
from torch_geometric.data import HeteroData


# ─────────────────────────────────────────────────────────────────────────────
# ENCODERS
# ─────────────────────────────────────────────────────────────────────────────

class SAGEEncoder(nn.Module):
    """
    GraphSAGE encoder for heterogeneous bipartite graph.
    Aggregates neighbour embeddings via mean pooling.

    Args:
        in_dim    : input feature dimension (320 for ESM2 8M)
        hidden_dim: hidden layer dimension
        out_dim   : output embedding dimension
        n_layers  : number of SAGE convolution layers
        dropout   : dropout rate between layers
    """
    def __init__(self, in_dim: int, hidden_dim: int, out_dim: int,
                 n_layers: int = 2, dropout: float = 0.1):
        super().__init__()
        self.dropout = dropout

        dims = [in_dim] + [hidden_dim] * (n_layers - 1) + [out_dim]
        self.convs = nn.ModuleList([
            SAGEConv(dims[i], dims[i + 1])
            for i in range(n_layers)
        ])
        self.norms = nn.ModuleList([
            nn.LayerNorm(dims[i + 1])
            for i in range(n_layers)
        ])

    def forward(self, x, edge_index):
        # x can be a tuple (x_src, x_dst) for bipartite graphs
        # SAGEConv handles this natively — output shape = x_dst.size(0)
        for i, (conv, norm) in enumerate(zip(self.convs, self.norms)):
            x_dst = x[1] if isinstance(x, tuple) else x
            x = conv(x, edge_index)
            x = norm(x)
            if i < len(self.convs) - 1:
                x = F.gelu(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
                # after first layer, x is now a plain tensor (dst nodes only)
        return x


class GATEncoder(nn.Module):
    """
    Graph Attention Network encoder for heterogeneous bipartite graph.
    Attention weights indicate which neighbours are most relevant.

    Args:
        in_dim    : input feature dimension
        hidden_dim: hidden layer dimension
        out_dim   : output embedding dimension
        n_layers  : number of GAT convolution layers
        n_heads   : number of attention heads per layer
        dropout   : dropout rate
    """
    def __init__(self, in_dim: int, hidden_dim: int, out_dim: int,
                 n_layers: int = 2, n_heads: int = 4, dropout: float = 0.1):
        super().__init__()
        self.dropout = dropout

        self.convs = nn.ModuleList()
        self.norms = nn.ModuleList()

        for i in range(n_layers):
            in_ch  = in_dim if i == 0 else hidden_dim * n_heads
            out_ch = hidden_dim if i < n_layers - 1 else out_dim
            heads  = n_heads  if i < n_layers - 1 else 1
            concat = True     if i < n_layers - 1 else False

            self.convs.append(GATConv(in_ch, out_ch, heads=heads,
                                       concat=concat, dropout=dropout,
                                       add_self_loops=False))
            self.norms.append(nn.LayerNorm(out_ch * heads if concat else out_ch))

    def forward(self, x, edge_index):
        # x can be a tuple (x_src, x_dst) for bipartite graphs
        # GATConv handles this natively — output shape = x_dst.size(0)
        for i, (conv, norm) in enumerate(zip(self.convs, self.norms)):
            x = conv(x, edge_index)
            x = norm(x)
            if i < len(self.convs) - 1:
                x = F.gelu(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
        return x


# ─────────────────────────────────────────────────────────────────────────────
# INPUT PROJECTION
# ─────────────────────────────────────────────────────────────────────────────

class InputProjection(nn.Module):
    """
    Project raw ESM2 features to a common hidden dimension before GNN layers.
    This allows TCR and epitope features (both 320-dim) to be projected
    independently before message passing.
    """
    def __init__(self, in_dim: int, out_dim: int, dropout: float = 0.1):
        super().__init__()
        self.proj = nn.Sequential(
            nn.Linear(in_dim, out_dim),
            nn.LayerNorm(out_dim),
            nn.GELU(),
            nn.Dropout(dropout),
        )

    def forward(self, x):
        return self.proj(x)


# ─────────────────────────────────────────────────────────────────────────────
# LINK PREDICTION HEADS
# ─────────────────────────────────────────────────────────────────────────────

class DotProductDecoder(nn.Module):
    """
    Score a (TCR, epitope) pair by dot product of their embeddings.
    Fast, parameter-free, symmetric.
    """
    def forward(self, tcr_emb, epi_emb):
        return (tcr_emb * epi_emb).sum(dim=-1)   # (E,)


class MLPDecoder(nn.Module):
    """
    Score a (TCR, epitope) pair via MLP on concatenated embeddings.
    More expressive than dot product; can learn asymmetric interactions.

    Args:
        in_dim  : 2 * node_embedding_dim
        hidden_dim: MLP hidden size
        dropout : dropout rate
    """
    def __init__(self, in_dim: int, hidden_dim: int = 128, dropout: float = 0.1):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1),
        )

    def forward(self, tcr_emb, epi_emb):
        x = torch.cat([tcr_emb, epi_emb], dim=-1)   # (E, 2*D)
        return self.mlp(x).squeeze(-1)               # (E,)


# ─────────────────────────────────────────────────────────────────────────────
# FULL MODEL
# ─────────────────────────────────────────────────────────────────────────────

class TCREpitopeGNN(nn.Module):
    """
    Full bipartite GNN for TCR–Epitope binding link prediction.

    Forward pass:
      1. Project TCR and epitope features independently
      2. Run GNN encoder (SAGE or GAT) on the bipartite graph
      3. Decode (TCR_i, Epitope_j) pairs with dot product or MLP
      4. Return logits (binary cross-entropy applied outside)

    Args:
        esm_dim     : ESM2 embedding dimension (320 for 8M model)
        hidden_dim  : GNN hidden dimension
        out_dim     : final node embedding dimension
        n_layers    : GNN depth
        encoder_type: "sage" or "gat"
        decoder_type: "dot" or "mlp"
        n_heads     : attention heads (GAT only)
        dropout     : dropout rate throughout
    """
    def __init__(
        self,
        esm_dim     : int   = 320,
        hidden_dim  : int   = 128,
        out_dim     : int   = 64,
        n_layers    : int   = 2,
        encoder_type: str   = "sage",
        decoder_type: str   = "dot",
        n_heads     : int   = 4,
        dropout     : float = 0.1,
    ):
        super().__init__()
        self.encoder_type = encoder_type
        self.decoder_type = decoder_type

        # independent input projections for TCR and epitope
        self.tcr_proj = InputProjection(esm_dim, hidden_dim, dropout)
        self.epi_proj = InputProjection(esm_dim, hidden_dim, dropout)

        # shared GNN encoder
        if encoder_type == "sage":
            self.encoder = SAGEEncoder(hidden_dim, hidden_dim, out_dim,
                                        n_layers, dropout)
        elif encoder_type == "gat":
            self.encoder = GATEncoder(hidden_dim, hidden_dim, out_dim,
                                       n_layers, n_heads, dropout)
        else:
            raise ValueError(f"encoder_type must be 'sage' or 'gat', got '{encoder_type}'")

        # link prediction decoder
        if decoder_type == "dot":
            self.decoder = DotProductDecoder()
        elif decoder_type == "mlp":
            self.decoder = MLPDecoder(out_dim * 2, hidden_dim, dropout)
        else:
            raise ValueError(f"decoder_type must be 'dot' or 'mlp', got '{decoder_type}'")

    def encode(self, data: HeteroData):
        """
        Run input projection + GNN encoder on the heterogeneous bipartite graph.
        Returns (tcr_embeddings, epitope_embeddings).

        For bipartite graphs, SAGEConv and GATConv accept a tuple (x_src, x_dst)
        where indices in edge_index refer to positions within x_src and x_dst
        respectively — never concatenate node tensors or offset indices.

        Two message-passing passes:
          1. TCR -> Epitope: epitopes aggregate from TCR neighbours
          2. Epitope -> TCR: TCRs aggregate from epitope neighbours
        """
        tcr_x = self.tcr_proj(data["tcr"].x)      # (N_tcr, hidden_dim)
        epi_x = self.epi_proj(data["epitope"].x)  # (N_epi, hidden_dim)

        edge_fwd = data["tcr", "binds",      "epitope"].edge_index  # src=tcr, dst=epi
        edge_rev = data["epitope", "rev_binds", "tcr"].edge_index   # src=epi, dst=tcr

        # epitopes aggregate from TCR neighbours
        epi_out = self.encoder((tcr_x, epi_x), edge_fwd)  # (N_epi, out_dim)

        # TCRs aggregate from epitope neighbours
        tcr_out = self.encoder((epi_x, tcr_x), edge_rev)  # (N_tcr, out_dim)

        return tcr_out, epi_out

    def decode(self, tcr_emb, epi_emb, edge_label_index):
        """
        Score candidate (TCR, Epitope) pairs.

        Args:
            tcr_emb         : (N_tcr, out_dim)
            epi_emb         : (N_epi, out_dim)
            edge_label_index: (2, E_candidates) — pairs to score

        Returns:
            logits: (E_candidates,)
        """
        src = edge_label_index[0]   # TCR indices
        dst = edge_label_index[1]   # epitope indices
        return self.decoder(tcr_emb[src], epi_emb[dst])

    def forward(self, data: HeteroData, edge_label_index):
        tcr_emb, epi_emb = self.encode(data)
        return self.decode(tcr_emb, epi_emb, edge_label_index)

    def count_parameters(self):
        return sum(p.numel() for p in self.parameters() if p.requires_grad)


# ─────────────────────────────────────────────────────────────────────────────
# FACTORY
# ─────────────────────────────────────────────────────────────────────────────

def build_model(cfg: dict) -> TCREpitopeGNN:
    """Build a TCREpitopeGNN from a config dict (loaded from YAML)."""
    model_cfg = cfg["model"]
    return TCREpitopeGNN(
        esm_dim      = model_cfg.get("esm_dim",      320),
        hidden_dim   = model_cfg.get("hidden_dim",   128),
        out_dim      = model_cfg.get("out_dim",       64),
        n_layers     = model_cfg.get("n_layers",       2),
        encoder_type = model_cfg.get("encoder_type", "sage"),
        decoder_type = model_cfg.get("decoder_type", "dot"),
        n_heads      = model_cfg.get("n_heads",        4),
        dropout      = model_cfg.get("dropout",       0.1),
    )