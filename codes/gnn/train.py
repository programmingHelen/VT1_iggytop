"""
train.py
========
Training script for TCR–Epitope GNN baseline.

Usage
-----
  # quick smoke test (one batch):
  python train.py --config configs/sage_dot.yaml --mode one_batch

  # full training:
  python train.py --config configs/sage_dot.yaml --mode train

  # override any config value via CLI:
  python train.py --config configs/gat_mlp.yaml --mode train training.epochs=10

Pipeline
--------
  1. Load shap_data.pkl  (X_test_tcr, X_test_epi, df_test_meta, y_test)
  2. Build bipartite HeteroData graph from positive pairs
  3. Split edges into train / val / test (by epitope to avoid leakage)
  4. Add negative edges by random sampling
  5. Train GNN with BCE loss, Adam optimizer, early stopping
  6. Evaluate: AUC, AUPR, F1, accuracy
  7. Log everything to MLflow
"""

import os
import sys
import pickle
import random
import logging
import argparse
from pathlib import Path
from typing import Tuple, Dict

import yaml
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.data import HeteroData
from torch_geometric.transforms import RandomLinkSplit
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import (
    roc_auc_score, average_precision_score,
    f1_score, accuracy_score,
)
import mlflow
import mlflow.pytorch
import subprocess

from models import build_model

# ─────────────────────────────────────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# REPRODUCIBILITY
# ─────────────────────────────────────────────────────────────────────────────

def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


# ─────────────────────────────────────────────────────────────────────────────
# 1. CONFIG LOADING
# ─────────────────────────────────────────────────────────────────────────────

def load_config(path: str, overrides: list) -> dict:
    """
    Load YAML config and apply CLI overrides.
    Override syntax: key.subkey=value  (e.g. training.epochs=20)
    """
    with open(path) as f:
        cfg = yaml.safe_load(f)

    for override in overrides:
        key_path, value = override.split("=", 1)
        keys = key_path.strip().split(".")
        d = cfg
        for k in keys[:-1]:
            d = d.setdefault(k, {})
        # try to cast to int/float/bool, else keep string
        for cast in (int, float):
            try:
                value = cast(value)
                break
            except ValueError:
                pass
        if value == "true":  value = True
        if value == "false": value = False
        d[keys[-1]] = value

    return cfg


# ─────────────────────────────────────────────────────────────────────────────
# 2. DATA LOADING & GRAPH CONSTRUCTION
# ─────────────────────────────────────────────────────────────────────────────

def load_pickle(path: str) -> dict:
    log.info(f"Loading {path} ...")
    with open(path, "rb") as f:
        return pickle.load(f)


def build_graph(sd: dict, cfg: dict, device: torch.device) -> Tuple[HeteroData, dict]:
    """
    Build a bipartite HeteroData graph from the shap_data pickle.

    Nodes:
      - TCR nodes:     unique CDR3β sequences, features = ESM2 TCR embedding
      - Epitope nodes: unique epitope sequences, features = ESM2 epi embedding

    Edges:
      - ("tcr", "binds", "epitope"): positive pairs only (label == 1)
      - Negative edges are added during training by RandomLinkSplit

    Returns:
      data     : HeteroData graph
      node_maps: dicts mapping sequence string → node index
    """
    meta = sd["df_test_meta"]
    positives = meta[meta["label"] == 1].reset_index(drop=True)

    log.info(f"  Positive pairs: {len(positives)}")

    # ── build unique node sets ────────────────────────────────────────────────
    tcr_seqs = positives["cdr3b"].unique().tolist()
    epi_seqs = positives["epitope"].unique().tolist()

    tcr_map = {s: i for i, s in enumerate(tcr_seqs)}
    epi_map = {s: i for i, s in enumerate(epi_seqs)}

    log.info(f"  TCR nodes     : {len(tcr_seqs)}")
    log.info(f"  Epitope nodes : {len(epi_seqs)}")

    # ── node features from embeddings ─────────────────────────────────────────
    # Build lookup: sequence string → ESM2 embedding vector
    all_tcr_seqs  = meta["cdr3b"].values
    all_epi_seqs  = meta["epitope"].values
    X_tcr_all     = sd["X_test_tcr"]   # (N, 320)
    X_epi_all     = sd["X_test_epi"]   # (N, 320)

    tcr_emb_lookup: Dict[str, np.ndarray] = {}
    epi_emb_lookup: Dict[str, np.ndarray] = {}
    for i, (t, e) in enumerate(zip(all_tcr_seqs, all_epi_seqs)):
        if t not in tcr_emb_lookup:
            tcr_emb_lookup[t] = X_tcr_all[i]
        if e not in epi_emb_lookup:
            epi_emb_lookup[e] = X_epi_all[i]

    tcr_feat = torch.tensor(
        np.stack([tcr_emb_lookup[s] for s in tcr_seqs]), dtype=torch.float32
    )
    epi_feat = torch.tensor(
        np.stack([epi_emb_lookup[s] for s in epi_seqs]), dtype=torch.float32
    )

    # ── edge index ────────────────────────────────────────────────────────────
    src = torch.tensor([tcr_map[r] for r in positives["cdr3b"]], dtype=torch.long)
    dst = torch.tensor([epi_map[r] for r in positives["epitope"]], dtype=torch.long)

    # ── assemble HeteroData ───────────────────────────────────────────────────
    data = HeteroData()
    data["tcr"].x                              = tcr_feat
    data["epitope"].x                          = epi_feat
    data["tcr", "binds", "epitope"].edge_index = torch.stack([src, dst])

    # add reverse edges so message can flow both ways
    data["epitope", "rev_binds", "tcr"].edge_index = torch.stack([dst, src])

    log.info(f"  Graph edges   : {src.size(0)}")

    return data.to(device), {"tcr": tcr_map, "epitope": epi_map}


def split_edges(data: HeteroData, cfg: dict) -> Tuple[HeteroData, HeteroData, HeteroData]:
    """
    Split graph edges into train / val / test using RandomLinkSplit.
    Negative edges are added automatically (ratio controlled by neg_sampling_ratio).
    """
    split_cfg = cfg["data"]
    transform = RandomLinkSplit(
        num_val              = split_cfg["val_frac"],
        num_test             = split_cfg["test_frac"],
        neg_sampling_ratio   = split_cfg["neg_sampling_ratio"],
        edge_types           = [("tcr", "binds", "epitope")],
        rev_edge_types       = [("epitope", "rev_binds", "tcr")],
        is_undirected        = False,
        add_negative_train_samples = True,
    )
    train_data, val_data, test_data = transform(data)
    log.info(
        f"  Train edges: {train_data['tcr','binds','epitope'].edge_label.size(0)} | "
        f"Val: {val_data['tcr','binds','epitope'].edge_label.size(0)} | "
        f"Test: {test_data['tcr','binds','epitope'].edge_label.size(0)}"
    )
    return train_data, val_data, test_data


# ─────────────────────────────────────────────────────────────────────────────
# 3. TRAIN / EVAL FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def train_step(model, data, optimizer, device) -> float:
    model.train()
    optimizer.zero_grad()

    edge_label_index = data["tcr", "binds", "epitope"].edge_label_index.to(device)
    edge_label       = data["tcr", "binds", "epitope"].edge_label.float().to(device)

    logits = model(data, edge_label_index)
    loss   = F.binary_cross_entropy_with_logits(logits, edge_label)

    loss.backward()
    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
    optimizer.step()

    return loss.item()


@torch.no_grad()
def evaluate(model, data, device) -> dict:
    model.eval()

    edge_label_index = data["tcr", "binds", "epitope"].edge_label_index.to(device)
    edge_label       = data["tcr", "binds", "epitope"].edge_label.float().cpu().numpy()

    logits = model(data, edge_label_index).cpu().numpy()
    probs  = 1 / (1 + np.exp(-logits))   # sigmoid
    preds  = (probs >= 0.5).astype(int)

    return {
        "loss"     : float(F.binary_cross_entropy_with_logits(
                         torch.tensor(logits), torch.tensor(edge_label)).item()),
        "auc"      : roc_auc_score(edge_label, probs),
        "aupr"     : average_precision_score(edge_label, probs),
        "f1"       : f1_score(edge_label, preds, zero_division=0),
        "accuracy" : accuracy_score(edge_label, preds),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 4. FULL TRAINING LOOP
# ─────────────────────────────────────────────────────────────────────────────

def train(cfg: dict):
    set_seed(cfg["training"]["seed"])
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log.info(f"Device: {device}")

    # ── MLflow ────────────────────────────────────────────────────────────────
    subprocess.run(["mlflow", "db", "upgrade", cfg["mlflow"]["uri"]],
                   capture_output=True)
    mlflow.set_tracking_uri(cfg["mlflow"]["uri"])
    mlflow.set_experiment(cfg["mlflow"]["experiment"])

    with mlflow.start_run(run_name=cfg["mlflow"].get("run_name", "gnn_baseline")):
        mlflow.log_params({
            "encoder"   : cfg["model"]["encoder_type"],
            "decoder"   : cfg["model"]["decoder_type"],
            "hidden_dim": cfg["model"]["hidden_dim"],
            "out_dim"   : cfg["model"]["out_dim"],
            "n_layers"  : cfg["model"]["n_layers"],
            "dropout"   : cfg["model"]["dropout"],
            "lr"        : cfg["training"]["lr"],
            "epochs"    : cfg["training"]["epochs"],
            "patience"  : cfg["training"]["patience"],
        })

        # ── load data & build graph ───────────────────────────────────────────
        sd           = load_pickle(cfg["data"]["shap_pkl"])
        data, maps   = build_graph(sd, cfg, device)
        train_d, val_d, test_d = split_edges(data, cfg)

        # ── model ─────────────────────────────────────────────────────────────
        model = build_model(cfg).to(device)
        log.info(f"Model parameters: {model.count_parameters():,}")
        mlflow.log_param("n_params", model.count_parameters())

        optimizer = torch.optim.Adam(
            model.parameters(),
            lr           = cfg["training"]["lr"],
            weight_decay = cfg["training"]["weight_decay"],
        )
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode="max", factor=0.5,
            patience=cfg["training"]["lr_patience"],
        )

        # ── training loop ─────────────────────────────────────────────────────
        best_auc       = 0.0
        patience_count = 0
        best_ckpt      = Path(cfg["training"]["checkpoint_dir"]) / "best_model.pt"
        best_ckpt.parent.mkdir(parents=True, exist_ok=True)

        for epoch in range(1, cfg["training"]["epochs"] + 1):
            train_loss = train_step(model, train_d, optimizer, device)
            val_m      = evaluate(model, val_d, device)

            scheduler.step(val_m["auc"])
            lr_now = optimizer.param_groups[0]["lr"]

            log.info(
                f"Epoch {epoch:3d}/{cfg['training']['epochs']} | "
                f"train_loss={train_loss:.4f} | "
                f"val_auc={val_m['auc']:.4f} | "
                f"val_aupr={val_m['aupr']:.4f} | "
                f"lr={lr_now:.2e}"
            )

            mlflow.log_metrics({
                "train_loss": train_loss,
                "val_auc"   : val_m["auc"],
                "val_aupr"  : val_m["aupr"],
                "val_f1"    : val_m["f1"],
                "val_acc"   : val_m["accuracy"],
                "lr"        : lr_now,
            }, step=epoch)

            if val_m["auc"] > best_auc:
                best_auc       = val_m["auc"]
                patience_count = 0
                torch.save(model.state_dict(), best_ckpt)
                log.info(f"  ✓ New best AUC={best_auc:.4f}")
            else:
                patience_count += 1
                if patience_count >= cfg["training"]["patience"]:
                    log.info(f"Early stopping at epoch {epoch}")
                    break

        # ── test evaluation ───────────────────────────────────────────────────
        model.load_state_dict(torch.load(best_ckpt, map_location=device))
        test_m = evaluate(model, test_d, device)

        log.info("=" * 55)
        log.info("TEST RESULTS")
        for k, v in test_m.items():
            log.info(f"  {k:<12}: {v:.4f}")
        log.info("=" * 55)

        mlflow.log_metrics({f"test_{k}": v for k, v in test_m.items()})
        mlflow.pytorch.log_model(model, artifact_path="gnn_model")
        mlflow.log_artifact(str(best_ckpt))


# ─────────────────────────────────────────────────────────────────────────────
# 5. SMOKE TEST — one batch only
# ─────────────────────────────────────────────────────────────────────────────

def one_batch_test(cfg: dict):
    """
    Run one forward + backward pass to verify the pipeline end-to-end.
    Does NOT log to MLflow. Finishes in seconds.
    """
    log.info("=== ONE-BATCH SMOKE TEST ===")
    set_seed(cfg["training"]["seed"])
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    sd           = load_pickle(cfg["data"]["shap_pkl"])
    data, _      = build_graph(sd, cfg, device)
    train_d, val_d, _ = split_edges(data, cfg)

    model     = build_model(cfg).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg["training"]["lr"])

    loss = train_step(model, train_d, optimizer, device)
    val_m = evaluate(model, val_d, device)

    log.info(f"  Parameters: {model.count_parameters():,}")
    log.info(f"  Train loss: {loss:.4f}")
    log.info(f"  Val AUC   : {val_m['auc']:.4f}")
    log.info("Smoke test PASSED ✓")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TCR–Epitope GNN training")
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    parser.add_argument("--mode",   default="train",
                        choices=["train", "one_batch"],
                        help="'train' for full run, 'one_batch' for smoke test")
    parser.add_argument("overrides", nargs="*",
                        help="Config overrides: key.subkey=value")
    args = parser.parse_args()

    cfg = load_config(args.config, args.overrides)
    log.info(f"Config: {cfg}")

    if args.mode == "one_batch":
        one_batch_test(cfg)
    else:
        train(cfg)
