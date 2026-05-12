"""
GNN Models for TCR-Epitope Binding Prediction

Implementations of GraphSAGE, GAT, and GCN using PyTorch Geometric.
All models inherit from BaseGNNModel for consistent interface.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import (
    SAGEConv, GATConv, GCNConv,
    BatchNorm
)
from typing import Optional, List, Tuple
from abc import ABC, abstractmethod


class BaseGNNModel(nn.Module, ABC):
    """Abstract base class for all GNN models."""
    
    def __init__(
        self,
        input_dim: int,
        hidden_dim: int,
        output_dim: int,
        num_layers: int,
        dropout: float = 0.3,
        batch_norm: bool = True,
        activation: str = "relu",
        **kwargs
    ):
        super().__init__()
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.output_dim = output_dim
        self.num_layers = num_layers
        self.dropout = dropout
        self.batch_norm = batch_norm
        self.activation_name = activation
        
        # Activation function
        if activation == "relu":
            self.activation = F.relu
        elif activation == "elu":
            self.activation = F.elu
        elif activation == "gelu":
            self.activation = F.gelu
        else:
            self.activation = F.relu
    
    @abstractmethod
    def forward(self, x, edge_index, batch=None):
        """Forward pass. Subclasses must implement."""
        pass
    
    def _build_bn_layers(self, dims: list) -> nn.ModuleList:
        """Build batch norm layers if enabled.
        
        Args:
            dims: List of layer output dimensions [d1, d2, ..., d_n]
        """
        if not self.batch_norm:
            return nn.ModuleList()
        
        bn_layers = nn.ModuleList()
        # Skip the first dimension (input), apply BN to each layer output
        for d in dims[1:]:
            bn_layers.append(BatchNorm(d))
        return bn_layers


class GraphSAGE(BaseGNNModel):
    """
    GraphSAGE: Graph Inductive Representation Learning
    
    References:
        Hamilton et al., "Inductive Representation Learning on Large Graphs", NIPS 2017
    """
    
    def __init__(
        self,
        input_dim: int,
        hidden_dim: int,
        output_dim: int,
        num_layers: int,
        aggregator_type: str = "mean",
        sample_sizes: Optional[List[int]] = None,
        dropout: float = 0.3,
        batch_norm: bool = True,
        activation: str = "relu",
        **kwargs
    ):
        super().__init__(
            input_dim, hidden_dim, output_dim, num_layers,
            dropout, batch_norm, activation
        )
        
        self.aggregator_type = aggregator_type
        self.sample_sizes = sample_sizes or [10] * num_layers
        
        # Build layers with correct dimensions
        dims = [input_dim] + [hidden_dim] * (num_layers - 1) + [output_dim]
        self.layers = nn.ModuleList()
        self.bn_layers = self._build_bn_layers(dims)
        
        for i in range(num_layers):
            self.layers.append(
                SAGEConv(
                    dims[i], 
                    dims[i + 1],
                    aggr=aggregator_type
                )
            )
    
    def forward(self, x: torch.Tensor, edge_index: torch.Tensor, batch=None) -> torch.Tensor:
        """
        Args:
            x: Node feature matrix of shape [num_nodes, input_dim]
            edge_index: Edge indices of shape [2, num_edges]
            batch: Batch assignment vector (for mini-batch processing)
        
        Returns:
            Node embeddings of shape [num_nodes, output_dim]
        """
        for i, layer in enumerate(self.layers):
            x = layer(x, edge_index)
            
            if self.batch_norm and i < len(self.bn_layers):
                x = self.bn_layers[i](x)
            
            if i < len(self.layers) - 1:
                x = self.activation(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
        
        return x


class GAT(BaseGNNModel):
    """
    Graph Attention Network (GAT)
    
    References:
        Veličković et al., "Graph Attention Networks", ICLR 2018
    """
    
    def __init__(
        self,
        input_dim: int,
        hidden_dim: int,
        output_dim: int,
        num_layers: int,
        num_heads: int = 8,
        concat_heads: bool = True,
        dropout: float = 0.3,
        attention_dropout: float = 0.1,
        batch_norm: bool = True,
        activation: str = "relu",
        residual: bool = False,
        **kwargs
    ):
        super().__init__(
            input_dim, hidden_dim, output_dim, num_layers,
            dropout, batch_norm, activation
        )
        
        self.num_heads = num_heads
        self.concat_heads = concat_heads
        self.attention_dropout = attention_dropout
        self.residual = residual
        
        # Pre-compute all output dimensions for batch norm
        # Account for head concatenation: GATConv outputs out_dim * num_heads when concat=True
        dims = [input_dim]
        for i in range(num_layers):
            if i < num_layers - 1:
                # Intermediate layers: output is hidden_dim * num_heads if concat_heads
                dims.append(hidden_dim * (num_heads if concat_heads else 1))
            else:
                # Last layer: output is output_dim (heads=1)
                dims.append(output_dim)
        
        # Build layers
        self.layers = nn.ModuleList()
        self.bn_layers = self._build_bn_layers(dims)
        
        # For residual connections, need matching dimensions
        self.residual_layers = nn.ModuleList() if residual else None
        
        for i in range(num_layers):
            if i == 0:
                in_dim = input_dim
                out_dim = hidden_dim
            elif i == num_layers - 1:
                in_dim = hidden_dim * (num_heads if concat_heads else 1)
                out_dim = output_dim
            else:
                in_dim = hidden_dim * (num_heads if concat_heads else 1)
                out_dim = hidden_dim
            
            self.layers.append(
                GATConv(
                    in_dim,
                    out_dim,
                    heads=num_heads if i < num_layers - 1 else 1,
                    concat=concat_heads if i < num_layers - 1 else True,
                    dropout=attention_dropout,
                    add_self_loops=True
                )
            )
            
            # Residual projection if dimensions don't match
            if residual and i > 0:
                self.residual_layers.append(
                    nn.Linear(in_dim, out_dim) if in_dim != out_dim else nn.Identity()
                )
    
    def forward(self, x: torch.Tensor, edge_index: torch.Tensor, batch=None) -> torch.Tensor:
        """
        Args:
            x: Node feature matrix of shape [num_nodes, input_dim]
            edge_index: Edge indices of shape [2, num_edges]
            batch: Batch assignment vector
        
        Returns:
            Node embeddings of shape [num_nodes, output_dim]
        """
        for i, layer in enumerate(self.layers):
            x_in = x
            x = layer(x, edge_index)
            
            # Residual connection
            if self.residual and i > 0:
                x = x + self.residual_layers[i - 1](x_in)
            
            if self.batch_norm and i < len(self.bn_layers):
                x = self.bn_layers[i](x)
            
            if i < len(self.layers) - 1:
                x = self.activation(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
        
        return x


class GCN(BaseGNNModel):
    """
    Graph Convolutional Network (GCN)
    
    References:
        Kipf & Welling, "Semi-Supervised Classification with Graph Convolutional Networks", ICLR 2017
    """
    
    def __init__(
        self,
        input_dim: int,
        hidden_dim: int,
        output_dim: int,
        num_layers: int,
        dropout: float = 0.3,
        batch_norm: bool = True,
        bias: bool = True,
        activation: str = "relu",
        cached: bool = False,
        **kwargs
    ):
        super().__init__(
            input_dim, hidden_dim, output_dim, num_layers,
            dropout, batch_norm, activation
        )
        
        self.cached = cached
        
        # Build layers with correct dimensions
        dims = [input_dim] + [hidden_dim] * (num_layers - 1) + [output_dim]
        self.layers = nn.ModuleList()
        self.bn_layers = self._build_bn_layers(dims)
        
        for i in range(num_layers):
            self.layers.append(
                GCNConv(
                    dims[i],
                    dims[i + 1],
                    bias=bias,
                    cached=cached,
                    add_self_loops=True
                )
            )
    
    def forward(self, x: torch.Tensor, edge_index: torch.Tensor, batch=None) -> torch.Tensor:
        """
        Args:
            x: Node feature matrix of shape [num_nodes, input_dim]
            edge_index: Edge indices of shape [2, num_edges]
            batch: Batch assignment vector
        
        Returns:
            Node embeddings of shape [num_nodes, output_dim]
        """
        for i, layer in enumerate(self.layers):
            x = layer(x, edge_index)
            
            if self.batch_norm and i < len(self.bn_layers):
                x = self.bn_layers[i](x)
            
            if i < len(self.layers) - 1:
                x = self.activation(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
        
        return x


class GNNLinkPredictor(nn.Module):
    """
    Link prediction head for binding prediction.
    Combines two node embeddings and predicts binding probability.
    """
    
    def __init__(self, embedding_dim: int, hidden_dim: int = 64):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(embedding_dim * 2, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden_dim // 2, 1),
        )
    
    def forward(self, z1: torch.Tensor, z2: torch.Tensor) -> torch.Tensor:
        """
        Predict binding score between two nodes.
        
        Args:
            z1: Embeddings of first node type (TCR) [batch_size, embedding_dim]
            z2: Embeddings of second node type (Epitope) [batch_size, embedding_dim]
        
        Returns:
            Binding scores [batch_size, 1]
        """
        combined = torch.cat([z1, z2], dim=1)
        return torch.sigmoid(self.mlp(combined))


def build_model(
    model_name: str,
    input_dim: int,
    hidden_dim: int,
    output_dim: int,
    num_layers: int,
    **kwargs
) -> BaseGNNModel:
    """
    Factory function to build GNN models.
    
    Args:
        model_name: One of "graphsage", "gat", "gcn"
        input_dim: Input feature dimension
        hidden_dim: Hidden layer dimension
        output_dim: Output embedding dimension
        num_layers: Number of layers
        **kwargs: Model-specific hyperparameters
    
    Returns:
        Instantiated GNN model
    """
    model_name = model_name.lower()
    
    if model_name == "graphsage":
        return GraphSAGE(
            input_dim, hidden_dim, output_dim, num_layers, **kwargs
        )
    elif model_name == "gat":
        return GAT(
            input_dim, hidden_dim, output_dim, num_layers, **kwargs
        )
    elif model_name == "gcn":
        return GCN(
            input_dim, hidden_dim, output_dim, num_layers, **kwargs
        )
    else:
        raise ValueError(f"Unknown model: {model_name}")
