"""
Training script for GNN TCR-Epitope Binding Prediction

Usage:
    python -m src.iggytop.gnn.train --config graphsage --epochs 100 --lr 0.001
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, Tuple, Optional
import json
from datetime import datetime

import yaml
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.optim.lr_scheduler import CosineAnnealingLR, LinearLR, StepLR
from sklearn.metrics import (
    roc_auc_score, roc_curve, precision_recall_curve, average_precision_score,
    f1_score, precision_score, recall_score, accuracy_score, confusion_matrix
)
import mlflow
import mlflow.pytorch

from .models import build_model, GNNLinkPredictor
from .data_loader import BipartiteGraphDataset, create_data_loaders

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class GNNTrainer:
    """Main trainer class for GNN models."""
    
    def __init__(
        self,
        model_config: Dict,
        train_config: Dict,
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
        output_dir: str = "outputs/gnn_models",
        use_mlflow: bool = True,
    ):
        self.model_config = model_config
        self.train_config = train_config
        self.device = torch.device(device)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.use_mlflow = use_mlflow
        
        logger.info(f"Using device: {self.device}")
        
        # Initialize MLflow if enabled
        if self.use_mlflow:
            mlflow.set_tracking_uri(train_config.get("mlflow_uri"))
            mlflow.set_experiment(train_config.get("mlflow_experiment", "gnn_binding"))
    
    def train_epoch(
        self,
        model: nn.Module,
        predictor: nn.Module,
        train_loader,
        optimizer: optim.Optimizer,
        criterion: nn.Module,
    ) -> Dict[str, float]:
        """Train for one epoch."""
        model.train()
        predictor.train()
        
        total_loss = 0
        all_preds = []
        all_labels = []
        
        for batch in train_loader:
            batch = batch.to(self.device)
            optimizer.zero_grad()
            
            # Forward pass
            node_embeddings = model(batch.x, batch.edge_index, batch.batch)
            
            # Link prediction
            tcr_idx = batch.tcr_nodes if hasattr(batch, 'tcr_nodes') else batch.y_tcr
            epi_idx = batch.epitope_nodes if hasattr(batch, 'epitope_nodes') else batch.y_epi
            
            tcr_emb = node_embeddings[tcr_idx]
            epi_emb = node_embeddings[epi_idx]
            
            logits = predictor(tcr_emb, epi_emb)
            labels = batch.y.unsqueeze(1).float()
            
            loss = criterion(logits, labels)
            loss.backward()
            
            # Gradient clipping
            torch.nn.utils.clip_grad_norm_(
                list(model.parameters()) + list(predictor.parameters()),
                max_norm=1.0
            )
            
            optimizer.step()
            
            total_loss += loss.item()
            all_preds.extend(logits.detach().cpu().numpy().flatten())
            all_labels.extend(labels.detach().cpu().numpy().flatten())
        
        # Compute metrics
        all_preds = np.array(all_preds)
        all_labels = np.array(all_labels)
        
        metrics = {
            "loss": total_loss / len(train_loader),
            "auc": roc_auc_score(all_labels, all_preds) if len(np.unique(all_labels)) > 1 else 0.5,
            "aupr": average_precision_score(all_labels, all_preds),
        }
        
        return metrics
    
    @torch.no_grad()
    def evaluate(
        self,
        model: nn.Module,
        predictor: nn.Module,
        val_loader,
        criterion: nn.Module,
    ) -> Dict[str, float]:
        """Evaluate on validation set."""
        model.eval()
        predictor.eval()
        
        total_loss = 0
        all_preds = []
        all_labels = []
        
        for batch in val_loader:
            batch = batch.to(self.device)
            
            node_embeddings = model(batch.x, batch.edge_index, batch.batch)
            
            tcr_idx = batch.tcr_nodes if hasattr(batch, 'tcr_nodes') else batch.y_tcr
            epi_idx = batch.epitope_nodes if hasattr(batch, 'epitope_nodes') else batch.y_epi
            
            tcr_emb = node_embeddings[tcr_idx]
            epi_emb = node_embeddings[epi_idx]
            
            logits = predictor(tcr_emb, epi_emb)
            labels = batch.y.unsqueeze(1).float()
            
            loss = criterion(logits, labels)
            total_loss += loss.item()
            
            all_preds.extend(logits.cpu().numpy().flatten())
            all_labels.extend(labels.cpu().numpy().flatten())
        
        all_preds = np.array(all_preds)
        all_labels = np.array(all_labels)
        
        # Compute optimal F1 threshold
        precision, recall, thresholds = precision_recall_curve(all_labels, all_preds)
        f1_scores = 2 * (precision * recall) / (precision + recall + 1e-8)
        best_threshold = thresholds[np.argmax(f1_scores)]
        
        metrics = {
            "loss": total_loss / len(val_loader),
            "auc": roc_auc_score(all_labels, all_preds) if len(np.unique(all_labels)) > 1 else 0.5,
            "aupr": average_precision_score(all_labels, all_preds),
            "f1": f1_score(all_labels, (all_preds > best_threshold).astype(int)),
            "precision": precision_score(all_labels, (all_preds > best_threshold).astype(int), zero_division=0),
            "recall": recall_score(all_labels, (all_preds > best_threshold).astype(int), zero_division=0),
            "accuracy": accuracy_score(all_labels, (all_preds > best_threshold).astype(int)),
            "best_threshold": best_threshold,
        }
        
        return metrics
    
    def train(
        self,
        train_loader,
        val_loader,
        test_loader=None,
        num_epochs: Optional[int] = None,
    ):
        """Full training loop with early stopping."""
        num_epochs = num_epochs or self.train_config["num_epochs"]
        learning_rate = self.train_config["learning_rate"]
        weight_decay = float(self.train_config.get("weight_decay", 1e-5))
        patience = self.train_config.get("early_stopping_patience", 10)
        
        # Build model
        # Extract kwargs, avoiding duplicate parameters
        model_kwargs = {k: v for k, v in self.model_config.items() 
                       if k not in ["name", "hidden_dim", "output_dim", "num_layers"]}
        
        model = build_model(
            self.model_config["name"],
            input_dim=640,  # Concatenated TCR (320) + Epitope (320) embeddings
            hidden_dim=self.model_config.get("hidden_dim", 128),
            output_dim=self.model_config.get("output_dim", 64),
            num_layers=self.model_config.get("num_layers", 3),
            **model_kwargs
        ).to(self.device)
        
        predictor = GNNLinkPredictor(
            embedding_dim=self.model_config.get("output_dim", 64),
            hidden_dim=128
        ).to(self.device)
        
        logger.info(f"Model: {model}")
        logger.info(f"Predictor: {predictor}")
        
        # Optimizer and scheduler
        optimizer = optim.Adam(
            list(model.parameters()) + list(predictor.parameters()),
            lr=learning_rate,
            weight_decay=weight_decay
        )
        
        scheduler_name = self.train_config.get("scheduler", "cosine")
        if scheduler_name == "cosine":
            scheduler = CosineAnnealingLR(optimizer, T_max=num_epochs)
        elif scheduler_name == "linear":
            scheduler = LinearLR(optimizer, start_factor=1.0, end_factor=0.1, total_iters=num_epochs)
        elif scheduler_name == "step":
            scheduler = StepLR(optimizer, step_size=30, gamma=0.1)
        else:
            scheduler = None
        
        criterion = nn.BCELoss()
        
        # Start MLflow run
        if self.use_mlflow:
            mlflow.start_run()
            mlflow.log_params({
                "model": self.model_config["name"],
                "num_epochs": num_epochs,
                "learning_rate": learning_rate,
                "weight_decay": weight_decay,
                **self.model_config
            })
        
        best_val_auc = 0
        patience_counter = 0
        history = {"train": [], "val": []}
        
        try:
            for epoch in range(num_epochs):
                train_metrics = self.train_epoch(model, predictor, train_loader, optimizer, criterion)
                val_metrics = self.evaluate(model, predictor, val_loader, criterion)
                
                if scheduler:
                    scheduler.step()
                
                history["train"].append(train_metrics)
                history["val"].append(val_metrics)
                
                logger.info(
                    f"Epoch {epoch+1}/{num_epochs} | "
                    f"Train Loss: {train_metrics['loss']:.4f} | "
                    f"Val Loss: {val_metrics['loss']:.4f} | "
                    f"Val AUC: {val_metrics['auc']:.4f} | "
                    f"Val F1: {val_metrics['f1']:.4f}"
                )
                
                if self.use_mlflow:
                    mlflow.log_metrics({
                        **{f"train_{k}": v for k, v in train_metrics.items()},
                        **{f"val_{k}": v for k, v in val_metrics.items()},
                    }, step=epoch)
                
                # Early stopping
                if val_metrics["auc"] > best_val_auc:
                    best_val_auc = val_metrics["auc"]
                    patience_counter = 0
                    
                    # Save checkpoint
                    checkpoint_path = self.output_dir / f"{self.model_config['name']}_best.pt"
                    torch.save({
                        "epoch": epoch,
                        "model_state": model.state_dict(),
                        "predictor_state": predictor.state_dict(),
                        "optimizer_state": optimizer.state_dict(),
                        "val_metrics": val_metrics,
                    }, checkpoint_path)
                    logger.info(f"Saved best model to {checkpoint_path}")
                else:
                    patience_counter += 1
                    if patience_counter >= patience:
                        logger.info(f"Early stopping at epoch {epoch+1}")
                        break
            
            # Load best model
            checkpoint = torch.load(self.output_dir / f"{self.model_config['name']}_best.pt")
            model.load_state_dict(checkpoint["model_state"])
            predictor.load_state_dict(checkpoint["predictor_state"])
            
            # Test evaluation if provided
            if test_loader:
                test_metrics = self.evaluate(model, predictor, test_loader, criterion)
                logger.info(f"Test Metrics: {test_metrics}")
                if self.use_mlflow:
                    mlflow.log_metrics({f"test_{k}": v for k, v in test_metrics.items()})
            
            # Save training history
            history_path = self.output_dir / f"{self.model_config['name']}_history.json"
            with open(history_path, "w") as f:
                # Convert numpy values to Python floats for JSON serialization
                history_json = {
                    "train": [{k: float(v) for k, v in epoch.items()} for epoch in history["train"]],
                    "val": [{k: float(v) if isinstance(v, (np.integer, np.floating)) else v for k, v in epoch.items()} for epoch in history["val"]]
                }
                json.dump(history_json, f, indent=2)
            
            logger.info(f"Training history saved to {history_path}")
            
        finally:
            if self.use_mlflow:
                mlflow.end_run()
        
        return model, predictor, history


def load_config(config_path: str) -> Dict:
    """Load YAML configuration file."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description="Train GNN for TCR-Epitope binding prediction")
    parser.add_argument(
        "--config",
        type=str,
        choices=["graphsage", "gat", "gcn"],
        default="graphsage",
        help="Model configuration to use"
    )
    parser.add_argument("--epochs", type=int, help="Number of epochs")
    parser.add_argument("--lr", type=float, help="Learning rate")
    parser.add_argument("--batch-size", type=int, help="Batch size")
    parser.add_argument("--device", type=str, default="cuda", help="Device: cuda or cpu")
    parser.add_argument("--output-dir", type=str, default="outputs/gnn_models", help="Output directory")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--no-mlflow", action="store_true", help="Disable MLflow tracking")
    
    args = parser.parse_args()
    
    # Set seed
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    
    # Load configs
    config_dir = Path(__file__).parent / "config"
    model_config = load_config(config_dir / f"{args.config}.yaml")
    train_config = load_config(config_dir / "train.yaml")
    
    # Override with command-line args
    if args.epochs:
        model_config["training"]["num_epochs"] = args.epochs
    if args.lr:
        model_config["training"]["learning_rate"] = args.lr
    if args.batch_size:
        model_config["training"]["batch_size"] = args.batch_size
    
    # Merge training configs
    model_config["training"].update(train_config.get("optimizer", {}))
    
    logger.info(f"Model config: {model_config}")
    logger.info(f"Train config: {train_config}")
    
    # Create data loaders
    logger.info("Creating data loaders...")
    train_loader, val_loader, test_loader = create_data_loaders(
        shap_pkl_path="outputs/outputs_v2/shap_data.pkl",
        batch_size=model_config["training"]["batch_size"],
        val_split=model_config["training"]["val_split"],
        seed=args.seed,
    )
    
    # Initialize trainer
    trainer = GNNTrainer(
        model_config=model_config["model"],
        train_config=model_config["training"],
        device=args.device,
        output_dir=args.output_dir,
        use_mlflow=not args.no_mlflow,
    )
    
    # Train
    logger.info("Starting training...")
    model, predictor, history = trainer.train(
        train_loader,
        val_loader,
        test_loader,
        num_epochs=model_config["training"]["num_epochs"],
    )
    
    logger.info("Training completed!")


if __name__ == "__main__":
    main()
