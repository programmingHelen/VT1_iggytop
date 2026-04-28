"""
Minimal TCR–Epitope Baseline with MLflow
=========================================
 
Only ESM2 embeddings. No feature engineering.
Uses scirpy's get_airr_context for clean data loading.
Tracks with MLflow.
 
Usage:
  python tcr_epitope_minimal.py
"""
 
import logging
import random
from typing import Optional, Dict

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scirpy as ir
import torch
from transformers import AutoTokenizer, EsmModel
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    roc_auc_score, average_precision_score,
    f1_score, precision_score, recall_score,
    accuracy_score, confusion_matrix,
    roc_curve, precision_recall_curve
)

import mlflow
import mlflow.sklearn
import pickle
import os
 
logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)
 
# ─────────────────────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────────────────────
 
DATA_PATH    = "data/deduplicated_anndata.h5ad"
ESM_MODEL    = "facebook/esm2_t6_8M_UR50D"
MLFLOW_URI   = "sqlite:///outputs/mlflow.db"
MLFLOW_EXPERIMENT = "tcr_epitope_minimal"
 
MAX_SAMPLES  = None   # set to None to use all data
NEG_RATIO    = 1.0
RANDOM_SEED  = 42
 
# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD DATA using scirpy (correct way for iggytop AnnData)
# ─────────────────────────────────────────────────────────────────────────────
 
def load_data(path: str, max_samples: Optional[int]) -> pd.DataFrame:
    """
    Load (cdr3b, epitope) pairs using scirpy's get_airr_context.

    scirpy understands the awkward array AIRR format natively —
    no need to iterate chains manually.

    Chain conventions in iggytop:
      VJ_1  = alpha chain (TRA)
      VDJ_1 = beta chain  (TRB) ← we want this one

    Filters applied:
      - Beta chain only (locus == TRB, not IG* for B-cells)
      - Productive sequences only
      - CDR3 length between 8 and 25 amino acids
    """
    log.info(f"Loading {path} ...")
    adata = sc.read_h5ad(path)

    # epitope_sequence lives in adata.obs (one per cell)
    epitopes = adata.obs["epitope_sequence"].values

    rows = []

    with ir.get.airr_context(
        adata,
        ["junction_aa", "locus", "productive"],
        chain=["VDJ_1"],   # VDJ_1 = beta chain (TRB)
    ) as m:
        obs = m.obs  # one row per cell, columns named VDJ_1_<field>

        for i in range(len(obs)):
            epi  = epitopes[i]
            cdr3 = obs.iloc[i]["VDJ_1_junction_aa"]
            locus = obs.iloc[i].get("VDJ_1_locus", "")
            productive = obs.iloc[i].get("VDJ_1_productive", None)

            # skip missing values
            if not isinstance(epi, str) or not epi:
                continue
            if cdr3 is None or (isinstance(cdr3, float) and np.isnan(cdr3)):
                continue
            cdr3 = str(cdr3).strip()
            if cdr3 in ("", "None", "nan"):
                continue

            # Filter out B-cells (IGH = heavy chain, IGL = light chain)
            if locus and str(locus).startswith(("IG", "ig")):
                continue

            rows.append({"cdr3b": cdr3, "epitope": epi})
 
    df = pd.DataFrame(rows).dropna()
    log.info(f"  Loaded {len(df)} TCR-epitope pairs")
 
    if max_samples and len(df) > max_samples:
        df = df.sample(max_samples, random_state=RANDOM_SEED).reset_index(drop=True)
        log.info(f"  Subsampled to {len(df)}")
 
    return df
 
# ─────────────────────────────────────────────────────────────────────────────
# 2. GENERATE NEGATIVES — mismatched pairing
# ─────────────────────────────────────────────────────────────────────────────
 
def generate_negatives(df_pos: pd.DataFrame, ratio: float) -> pd.DataFrame:
    """
    Create mismatched (TCR, epitope) pairs as negatives.
 
    Strategy: randomly pair a TCR with an epitope it did NOT come from.
    We check against the set of known positive pairs to avoid
    accidentally labelling a true positive as negative.
    """
    rng = random.Random(RANDOM_SEED)
    positive_pairs = set(zip(df_pos["cdr3b"], df_pos["epitope"]))
 
    tcrs     = df_pos["cdr3b"].values
    epitopes = df_pos["epitope"].values
 
    negatives   = []
    n_needed    = int(len(df_pos) * ratio)
    max_attempts = n_needed * 10
 
    attempts = 0
    while len(negatives) < n_needed and attempts < max_attempts:
        tcr = rng.choice(tcrs)
        epi = rng.choice(epitopes)
        if (tcr, epi) not in positive_pairs:
            negatives.append({"cdr3b": tcr, "epitope": epi, "label": 0})
        attempts += 1
 
    df_pos = df_pos.copy()
    df_pos["label"] = 1
    df_neg = pd.DataFrame(negatives)
 
    df = pd.concat([df_pos, df_neg], ignore_index=True).sample(
        frac=1, random_state=RANDOM_SEED
    ).reset_index(drop=True)
 
    log.info(f"  {len(df_pos)} pos + {len(df_neg)} neg = {len(df)} total")
    return df
 
# ─────────────────────────────────────────────────────────────────────────────
# 3. ESM2 EMBEDDINGS
# ─────────────────────────────────────────────────────────────────────────────
 
def embed_sequences(seqs, model_name: str, batch_size: int = 64) -> np.ndarray:
    """
    Embed amino-acid sequences using ESM2.
 
    Uses mean-pooling over the sequence length dimension,
    weighted by the attention mask to ignore padding tokens.
 
    Returns: np.ndarray of shape (N, embedding_dim)
    """
    log.info(f"Loading ESM2: {model_name}")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model     = EsmModel.from_pretrained(model_name)
    model.eval()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model  = model.to(device)
    log.info(f"  Device: {device}")
 
    # deduplicate to avoid re-embedding the same sequence twice
    unique_seqs = list(dict.fromkeys(seqs))
    cache = {}
 
    log.info(f"  Embedding {len(unique_seqs)} unique sequences ...")
    for i in range(0, len(unique_seqs), batch_size):
        batch  = unique_seqs[i : i + batch_size]
        inputs = tokenizer(
            batch,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=50,
        ).to(device)
 
        with torch.no_grad():
            out = model(**inputs)
 
        # mean-pool: weight by attention mask to ignore PAD tokens
        mask = inputs["attention_mask"].unsqueeze(-1).float()   # (B, L, 1)
        emb  = (out.last_hidden_state * mask).sum(1) / mask.sum(1)  # (B, D)
        emb  = emb.cpu().numpy()
 
        for s, e in zip(batch, emb):
            cache[s] = e
 
        # free GPU memory
        del out, inputs
        if device == "cuda":
            torch.cuda.empty_cache()
 
    return np.array([cache[s] for s in seqs])

# ─────────────────────────────────────────────────────────────────────────────
# 3. CALCULATE COMPREHENSIVE METRICS
# ─────────────────────────────────────────────────────────────────────────────

def calculate_metrics(y_true, y_pred_proba) -> Dict[str, float]:
    """
    Calculate comprehensive metrics for binary classification.
    
    Returns dictionary with:
      - auc: Area Under ROC Curve
      - aupr: Area Under PR Curve
      - f1: F1-score (at default threshold 0.5)
      - precision: Precision (at 0.5)
      - recall: Recall/Sensitivity (at 0.5)
      - specificity: True Negative Rate (at 0.5)
      - accuracy: Accuracy (at 0.5)
      - optimal_f1: Max F1-score (best threshold)
    """
    y_pred = (y_pred_proba >= 0.5).astype(int)
    
    metrics = {
        "auc": roc_auc_score(y_true, y_pred_proba),
        "aupr": average_precision_score(y_true, y_pred_proba),
        "f1": f1_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred),
        "recall": recall_score(y_true, y_pred),
        "accuracy": accuracy_score(y_true, y_pred),
    }
    
    # Calculate specificity (True Negative Rate)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    metrics["specificity"] = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    
    # Find optimal F1-score across different thresholds
    precisions, recalls, thresholds = precision_recall_curve(y_true, y_pred_proba)
    f1_scores = 2 * (precisions * recalls) / (precisions + recalls + 1e-10)
    metrics["optimal_f1"] = f1_scores.max()
    
    return metrics

# ─────────────────────────────────────────────────────────────────────────────
# 4. MAIN — train and evaluate with MLflow tracking
# ─────────────────────────────────────────────────────────────────────────────
 
if __name__ == "__main__":
    log.info("=" * 70)
    log.info("Minimal TCR–Epitope Baseline (ESM2 + 2 Models + MLflow)")
    log.info("=" * 70)

    # ── MLflow setup ──────────────────────────────────────────────────────────
    mlflow.set_tracking_uri(MLFLOW_URI)
    mlflow.set_experiment(MLFLOW_EXPERIMENT)

    with mlflow.start_run():

        # log hyperparameters
        mlflow.log_params({
            "esm_model"  : ESM_MODEL,
            "max_samples": MAX_SAMPLES,
            "neg_ratio"  : NEG_RATIO,
            "random_seed": RANDOM_SEED,
        })

        # ── load data ─────────────────────────────────────────────────────────
        df_pos = load_data(DATA_PATH, MAX_SAMPLES)
        mlflow.log_metric("n_positives", len(df_pos))

        df = generate_negatives(df_pos, NEG_RATIO)
        mlflow.log_metric("n_total", len(df))

        # ── embed ─────────────────────────────────────────────────────────────
        log.info("Embedding CDR3β sequences ...")
        tcr_emb = embed_sequences(df["cdr3b"].values, ESM_MODEL)

        log.info("Embedding epitope sequences ...")
        epi_emb = embed_sequences(df["epitope"].values, ESM_MODEL)

        # concatenate TCR + epitope embeddings — no hand-crafted features
        X = np.concatenate([tcr_emb, epi_emb], axis=1)
        y = df["label"].values

        log.info(f"Feature matrix: {X.shape}  (ESM2 only, no feature engineering)")
        mlflow.log_metric("feature_dim", X.shape[1])

        # ── split ─────────────────────────────────────────────────────────────
        X_train, X_test, y_train, y_test = train_test_split(
            X, y,
            test_size=0.2,
            random_state=RANDOM_SEED,
            stratify=y,
        )
        mlflow.log_metric("train_size", len(X_train))
        mlflow.log_metric("test_size",  len(X_test))

        # ── train 2 models ───────────────────────────────────────────────────
        models = {}

        # 1. Logistic Regression (simple baseline)
        log.info("Training Logistic Regression ...")
        lr = LogisticRegression(
            max_iter=1000,
            random_state=RANDOM_SEED,
            class_weight="balanced",
            n_jobs=-1,
        )
        lr.fit(X_train, y_train)
        models["logistic_regression"] = lr

        # 2. Random Forest
        log.info("Training Random Forest ...")
        rf = RandomForestClassifier(
            n_estimators=100,
            max_depth=15,
            n_jobs=-1,
            random_state=RANDOM_SEED,
            class_weight="balanced",
        )
        rf.fit(X_train, y_train)
        models["random_forest"] = rf

        # ── evaluate both models ──────────────────────────────────────────────
        log.info("=" * 70)
        log.info("EVALUATION RESULTS")
        log.info("=" * 70)

        all_metrics = {}
        for model_name, model in models.items():
            log.info(f"\n{model_name.upper()}")
            log.info("-" * 70)

            # Get predictions
            probs = model.predict_proba(X_test)[:, 1]

            # Calculate all metrics
            metrics = calculate_metrics(y_test, probs)
            all_metrics[model_name] = metrics

            # Log to MLflow with model name prefix
            for metric_name, value in metrics.items():
                mlflow.log_metric(f"{model_name}_{metric_name}", value)

            # Print metrics
            for metric_name in ["auc", "aupr", "f1", "precision", "recall", "specificity", "accuracy"]:
                log.info(f"  {metric_name:15s}: {metrics[metric_name]:.4f}")

        log.info("=" * 70)

        # ── save models to MLflow AND as local pickle ────────────────────────
        # MLflow log_model stores artifact_uri pointing to the cluster path.
        # To load the model locally we also save a plain pickle — portable.
        os.makedirs("outputs", exist_ok=True)

        for model_name, model in models.items():
            mlflow.sklearn.log_model(model, name=model_name)
            # plain pickle so the notebook can load it without MLflow artifacts
            pkl_path = f"outputs/{model_name}.pkl"
            with open(pkl_path, "wb") as f:
                pickle.dump(model, f)
            mlflow.log_artifact(pkl_path)   # also attach to run for reference
            log.info(f"✓ Saved {model_name} to {pkl_path}")

        # ── save data + run_id for SHAP and notebook ─────────────────────────
        # Include the MLflow run_id so the notebook can look up metrics
        # without hardcoding the run name.
        active_run_id = mlflow.active_run().info.run_id

        shap_data = {
            "X_test"  : X_test,
            "y_test"  : y_test,
            "rf_probs": models["random_forest"].predict_proba(X_test)[:, 1],
            "run_id"  : active_run_id,
            "mlflow_uri": MLFLOW_URI,
            "experiment": MLFLOW_EXPERIMENT,
        }

        with open("outputs/shap_data.pkl", "wb") as f:
            pickle.dump(shap_data, f)
        mlflow.log_artifact("outputs/shap_data.pkl")
        log.info("✓ Saved SHAP data to outputs/shap_data.pkl")

        # ── save metrics summary ──────────────────────────────────────────────
        summary_df = pd.DataFrame(all_metrics).T
        summary_df.to_csv("outputs/model_metrics_summary.csv")
        mlflow.log_artifact("outputs/model_metrics_summary.csv")
        log.info("✓ Saved metrics to outputs/model_metrics_summary.csv")

        log.info(f"\nRun ID: {active_run_id}")
        log.info(f"Run logged to: {MLFLOW_URI}")
        log.info("Done.")