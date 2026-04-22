"""
tcr_epitope_baseline.py
=======================
Classical ML baseline pipeline for TCR–epitope binding prediction.

Feature Engineering
-------------------
  - ESM2 embeddings of CDR3β and epitope sequences (facebook/esm2_t6_8M_UR50D)
  - Amino acid composition (AAC) for both sequences
  - Sequence lengths
  - V-gene one-hot encoding
  - 2-mer and 3-mer frequencies (optional, flag-controlled)

Negative Generation
-------------------
  Mismatched pairing strategy:
    For each positive (TCR_i, epitope_i), find a TCR from a *different* epitope
    group and pair it with epitope_i. This is biologically motivated:
    real non-binding pairs come from TCRs that bind other targets,
    not random shuffles which can accidentally create true positives.

Models
------
  - Logistic Regression
  - Random Forest

Experiment tracking: MLflow (auto-logged parameters, metrics, artifacts)

Usage
-----
  # One batch smoke test:
  python tcr_epitope_baseline.py --mode one_batch

  # Full run:
  python tcr_epitope_baseline.py --mode train

  # Skip ESM (use only hand-crafted features, faster):
  python tcr_epitope_baseline.py --mode train --no_esm
"""

# ─────────────────────────────────────────────────────────────────────────────
# IMPORTS
# ─────────────────────────────────────────────────────────────────────────────

import os
import gc
import math
import random
import logging
import argparse
import warnings
from pathlib import Path
from itertools import product
from collections import Counter
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import anndata as ad
import torch
from transformers import AutoTokenizer, EsmModel
import matplotlib.pyplot as plt

import mlflow
import mlflow.sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupShuffleSplit
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import (
    roc_auc_score, average_precision_score,
    accuracy_score, precision_score, recall_score,
    precision_recall_curve, roc_curve,
    confusion_matrix,
    classification_report,
)

warnings.filterwarnings("ignore", category=UserWarning)
# Suppress AnnData awkward array experimental warning (cosmetic, not a bug)
try:
    from anndata.utils import ExperimentalFeatureWarning
    warnings.filterwarnings("ignore", category=ExperimentalFeatureWarning)
except ImportError:
    pass

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────────────────────

DEFAULTS = dict(
    data_path    = "data/deduplicated_anndata.h5ad",
    output_dir   = "outputs",
    mlflow_uri   = "sqlite:///outputs/mlflow.db",  # SQLite backend (avoids FutureWarning)
    experiment   = "tcr_epitope_baselines",

    max_samples  = 10_000,    # cap positives loaded from AnnData (None = all)
    neg_ratio    = 1.0,       # negatives per positive
    use_kmers    = True,      # include 2-mer + 3-mer frequencies
    no_esm       = False,     # set True to skip ESM (much faster, weaker)
    esm_model    = "facebook/esm2_t6_8M_UR50D",
    esm_batch    = 64,        # sequences per ESM forward pass

    val_frac     = 0.15,
    test_frac    = 0.10,
    seed         = 42,
    n_seeds      = 1,       # run repeated experiments: seed, seed+1, ..., seed+n_seeds-1
    mode         = "train",   # "train" | "one_batch"

    # Logistic Regression
    lr_C         = 1.0,
    lr_max_iter  = 1000,

    # Random Forest
    rf_n_trees   = 200,
    rf_max_depth = 20,
    rf_n_jobs    = -1,
)

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

# ─────────────────────────────────────────────────────────────────────────────
# 1. DATA LOADING
# ─────────────────────────────────────────────────────────────────────────────

def load_from_anndata(path: str, max_samples: Optional[int]) -> pd.DataFrame:
    """
    Extract (cdr3b, epitope, v_gene) triples from AnnData.

    The AIRR field contains awkward.Record objects — NOT plain dicts.
    We use chain.tolist() to convert each record to a plain Python dict,
    then access fields with direct bracket notation (no .get()).

    Available chain fields confirmed by inspection:
      locus, productive, junction_aa, v_call, j_call, d_call,
      sequence, sequence_id, junction, ...
    """
    log.info(f"Loading AnnData from {path} …")
    adata    = ad.read_h5ad(path)
    airr     = adata.obsm["airr"]
    epitopes = adata.obs["epitope_sequence"]

    rows = []
    for i, chains in enumerate(airr):
        epi = epitopes.iloc[i]
        if not isinstance(epi, str) or not epi:
            continue
        for chain in chains:
            if chain is None:
                continue
            try:
                # Convert awkward Record → plain Python dict (the only safe way)
                c = chain.tolist()

                if c["locus"] != "TRB":
                    continue
                if not c["productive"]:
                    continue

                cdr3   = c["junction_aa"]
                v_call = c["v_call"]

                if not cdr3 or not v_call:
                    continue

                v_gene = v_call.split("*")[0].strip()

                if 8 <= len(cdr3) <= 25 and v_gene.startswith("TRBV"):
                    rows.append({"cdr3b": cdr3, "epitope": epi, "v_gene": v_gene})

            except Exception:
                continue

    df = pd.DataFrame(rows).dropna()
    log.info(f"Extracted {len(df)} positive (cdr3b, epitope, v_gene) rows")

    if max_samples and len(df) > max_samples:
        df = df.sample(max_samples, random_state=42).reset_index(drop=True)
        log.info(f"Subsampled to {len(df)}")

    return df


# ─────────────────────────────────────────────────────────────────────────────
# 2. NEGATIVE GENERATION — mismatched pairing
# ─────────────────────────────────────────────────────────────────────────────

def generate_negatives(df_pos: pd.DataFrame, neg_ratio: float, seed: int) -> pd.DataFrame:
    """
    Mismatched-pairing negatives:
      For each positive (TCR_i, epitope_j), sample a TCR from a *different*
      epitope group and pair it with epitope_j.

    This is better than random shuffling because:
      - Shuffling can accidentally re-create true positives
      - Mismatched pairing respects epitope identity:
        the negative TCR has known binding specificity to a different target

    Returns a DataFrame with the same columns + 'label' (0 for neg, 1 for pos).
    """
    rng = random.Random(seed)

    # group TCRs by epitope
    epitope_to_tcrs: dict = df_pos.groupby("epitope")["cdr3b"].apply(list).to_dict()
    epitope_to_vgene: dict = df_pos.groupby("epitope")["v_gene"].apply(list).to_dict()
    all_epitopes = list(epitope_to_tcrs.keys())

    n_neg = int(len(df_pos) * neg_ratio)
    neg_rows = []

    epitope_list = df_pos["epitope"].tolist()

    for _ in range(n_neg):
        # pick a random epitope as the "anchor"
        anchor_epi = rng.choice(epitope_list)

        # pick a different epitope to donate a TCR
        donor_epi = rng.choice(all_epitopes)
        attempts  = 0
        while donor_epi == anchor_epi and attempts < 20:
            donor_epi = rng.choice(all_epitopes)
            attempts += 1

        # sample a TCR from the donor epitope group
        donor_tcrs   = epitope_to_tcrs[donor_epi]
        donor_vgenes = epitope_to_vgene[donor_epi]
        idx          = rng.randrange(len(donor_tcrs))

        neg_rows.append({
            "cdr3b"  : donor_tcrs[idx],
            "epitope": anchor_epi,
            "v_gene" : donor_vgenes[idx],
        })

    df_neg = pd.DataFrame(neg_rows)
    df_neg["label"] = 0
    df_pos          = df_pos.copy()
    df_pos["label"] = 1

    df_all = pd.concat([df_pos, df_neg], ignore_index=True).sample(
        frac=1, random_state=seed
    ).reset_index(drop=True)

    log.info(
        f"Dataset: {len(df_pos)} positives + {len(df_neg)} negatives = {len(df_all)} total"
    )
    return df_all


# ─────────────────────────────────────────────────────────────────────────────
# 3. FEATURE ENGINEERING
# ─────────────────────────────────────────────────────────────────────────────

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

def amino_acid_composition(seq: str) -> np.ndarray:
    """20-dim vector: frequency of each standard amino acid."""
    counts = Counter(seq)
    total  = max(len(seq), 1)
    return np.array([counts.get(aa, 0) / total for aa in AMINO_ACIDS], dtype=np.float32)


def kmer_frequencies(seq: str, k: int) -> np.ndarray:
    """Frequency of all k-mers over the standard AA alphabet."""
    all_kmers = ["".join(p) for p in product(AMINO_ACIDS, repeat=k)]
    kmer_idx  = {km: i for i, km in enumerate(all_kmers)}
    counts    = np.zeros(len(all_kmers), dtype=np.float32)
    n_kmers   = max(len(seq) - k + 1, 1)
    for i in range(len(seq) - k + 1):
        km = seq[i : i + k]
        if km in kmer_idx:
            counts[kmer_idx[km]] += 1
    return counts / n_kmers


def encode_vgene_onehot(v_genes: pd.Series) -> Tuple[np.ndarray, List[str]]:
    """One-hot encode TRBV gene calls. Returns (matrix, category_list)."""
    cats    = sorted(v_genes.unique())
    cat_idx = {c: i for i, c in enumerate(cats)}
    mat     = np.zeros((len(v_genes), len(cats)), dtype=np.float32)
    for i, v in enumerate(v_genes):
        if v in cat_idx:
            mat[i, cat_idx[v]] = 1.0
    return mat, cats


# ── ESM2 embedding ───────────────────────────────────────────────────────────

class ESMEmbedder:
    """
    Wraps ESM2 to embed amino-acid sequences in batches.
    Returns mean-pooled last hidden states as numpy arrays.
    """
    def __init__(self, model_name: str, batch_size: int = 64):
        log.info(f"Loading ESM2 model: {model_name}")
        self.tokenizer  = AutoTokenizer.from_pretrained(model_name)
        self.model      = EsmModel.from_pretrained(model_name)
        self.model.eval()
        self.device     = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model      = self.model.to(self.device)
        self.batch_size = batch_size
        self._cache: dict = {}
        log.info(f"ESM2 on device: {self.device}")

    def embed(self, seqs: List[str]) -> np.ndarray:
        """Embed a list of sequences; returns (N, D) numpy array."""
        all_embs = []
        uncached = [s for s in seqs if s not in self._cache]

        # embed uncached sequences in batches
        for i in range(0, len(uncached), self.batch_size):
            batch = uncached[i : i + self.batch_size]
            inputs = self.tokenizer(
                batch,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=64,
            ).to(self.device)

            with torch.no_grad():
                outputs = self.model(**inputs)

            # mean-pool over sequence length (ignore padding via attention mask)
            mask = inputs["attention_mask"].unsqueeze(-1).float()  # (B, L, 1)
            emb  = (outputs.last_hidden_state * mask).sum(1) / mask.sum(1)
            emb  = emb.cpu().numpy()

            for seq, e in zip(batch, emb):
                self._cache[seq] = e

            # free GPU memory between batches
            del outputs, inputs
            if self.device.type == "cuda":
                torch.cuda.empty_cache()

        all_embs = np.array([self._cache[s] for s in seqs])
        return all_embs


# ── Batched feature extraction ───────────────────────────────────────────────

def extract_features_batched(
    df         : pd.DataFrame,
    embedder   : Optional[ESMEmbedder],
    use_kmers  : bool,
    batch_size : int = 1024,
    vgene_cats : Optional[List[str]] = None,
) -> Tuple[np.ndarray, List[str]]:
    """
    Extract features in row-batches to avoid memory spikes.

    Features (concatenated):
      - AAC of cdr3b                   (20 dims)
      - AAC of epitope                 (20 dims)
      - length of cdr3b                (1 dim)
      - length of epitope              (1 dim)
      - v_gene one-hot                 (n_vgenes dims)
      - 2-mer freqs cdr3b + epitope    (optional, 400*2 dims)
      - 3-mer freqs cdr3b + epitope    (optional, 8000*2 dims — expensive!)
      - ESM2 embedding cdr3b           (optional, 320 dims for 8M model)
      - ESM2 embedding epitope         (optional, 320 dims)

    Returns (X, feature_names).
    """
    # one-hot encode v_gene across the full dataframe so columns are consistent
    if vgene_cats is None:
        vgene_mat, vgene_cats = encode_vgene_onehot(df["v_gene"])
    else:
        cat_idx  = {c: i for i, c in enumerate(vgene_cats)}
        vgene_mat = np.zeros((len(df), len(vgene_cats)), dtype=np.float32)
        for i, v in enumerate(df["v_gene"]):
            if v in cat_idx:
                vgene_mat[i, cat_idx[v]] = 1.0

    # ESM embeddings for the full sequences upfront (already batched internally)
    if embedder is not None:
        log.info("Computing ESM2 embeddings for CDR3β sequences …")
        tcr_esm = embedder.embed(df["cdr3b"].tolist())
        log.info("Computing ESM2 embeddings for epitope sequences …")
        epi_esm = embedder.embed(df["epitope"].tolist())
    else:
        tcr_esm = epi_esm = None

    feature_blocks = []
    feat_names     = []

    n = len(df)
    for start in range(0, n, batch_size):
        end   = min(start + batch_size, n)
        batch = df.iloc[start:end]

        # amino acid composition
        tcr_aac = np.vstack([amino_acid_composition(s) for s in batch["cdr3b"]])
        epi_aac = np.vstack([amino_acid_composition(s) for s in batch["epitope"]])

        # lengths (normalised by typical max)
        tcr_len = (batch["cdr3b"].str.len().values / 25.0).reshape(-1, 1).astype(np.float32)
        epi_len = (batch["epitope"].str.len().values / 15.0).reshape(-1, 1).astype(np.float32)

        # v_gene one-hot slice
        vg_slice = vgene_mat[start:end]

        parts = [tcr_aac, epi_aac, tcr_len, epi_len, vg_slice]

        if use_kmers:
            # 2-mers only (3-mers produce 8000 dims per sequence — skip by default)
            tcr_2mer = np.vstack([kmer_frequencies(s, 2) for s in batch["cdr3b"]])
            epi_2mer = np.vstack([kmer_frequencies(s, 2) for s in batch["epitope"]])
            parts += [tcr_2mer, epi_2mer]

        if tcr_esm is not None:
            parts += [tcr_esm[start:end], epi_esm[start:end]]

        feature_blocks.append(np.concatenate(parts, axis=1))

    X = np.vstack(feature_blocks)

    # build feature names (only once)
    feat_names += [f"tcr_aac_{aa}" for aa in AMINO_ACIDS]
    feat_names += [f"epi_aac_{aa}" for aa in AMINO_ACIDS]
    feat_names += ["tcr_len", "epi_len"]
    feat_names += [f"vgene_{c}" for c in vgene_cats]
    if use_kmers:
        kmers = ["".join(p) for p in product(AMINO_ACIDS, repeat=2)]
        feat_names += [f"tcr_2mer_{k}" for k in kmers]
        feat_names += [f"epi_2mer_{k}" for k in kmers]
    if tcr_esm is not None:
        d = tcr_esm.shape[1]
        feat_names += [f"tcr_esm_{i}" for i in range(d)]
        feat_names += [f"epi_esm_{i}" for i in range(d)]

    log.info(f"Feature matrix shape: {X.shape}  ({len(feat_names)} features)")
    return X, feat_names, vgene_cats


# ─────────────────────────────────────────────────────────────────────────────
# 4. TRAIN / TEST SPLIT  (by epitope to avoid leakage)
# ─────────────────────────────────────────────────────────────────────────────

def split_by_epitope(
    df     : pd.DataFrame,
    X      : np.ndarray,
    test_frac : float,
    val_frac  : float,
    seed   : int,
) -> Tuple:
    """
    Split so that no epitope appears in both train and test.
    This is the correct way to evaluate generalisation to unseen epitopes.

    Limitation documented: epitope-based split may create class imbalance
    if some epitopes have many more samples than others.
    """
    groups = df["epitope"].values
    y      = df["label"].values

    gss_test = GroupShuffleSplit(n_splits=1, test_size=test_frac, random_state=seed)
    tv_idx, test_idx = next(gss_test.split(X, y, groups=groups))

    gss_val = GroupShuffleSplit(
        n_splits=1,
        test_size=val_frac / (1 - test_frac),
        random_state=seed,
    )
    train_idx, val_idx = next(
        gss_val.split(X[tv_idx], y[tv_idx], groups=groups[tv_idx])
    )
    train_idx = tv_idx[train_idx]
    val_idx   = tv_idx[val_idx]

    n_total = len(train_idx) + len(val_idx) + len(test_idx)
    log.info(
        f"Split (by epitope) → "
        f"train: {len(train_idx)} ({len(train_idx)/n_total:.1%}), "
        f"val: {len(val_idx)} ({len(val_idx)/n_total:.1%}), "
        f"test: {len(test_idx)} ({len(test_idx)/n_total:.1%})"
    )
    if len(test_idx) < 200:
        log.warning(
            f"Test set only has {len(test_idx)} samples — metrics will be noisy. "
            "Consider increasing max_samples or using a random split for baselines."
        )
    log.info(
        "Note: epitope-based split prevents leakage but may yield "
        "harder generalisation than random split."
    )
    return train_idx, val_idx, test_idx


# ─────────────────────────────────────────────────────────────────────────────
# 5. EVALUATION
# ─────────────────────────────────────────────────────────────────────────────

def choose_threshold_by_f1(y_true: np.ndarray, proba: np.ndarray) -> Tuple[float, float]:
    """Pick decision threshold that maximises F1 on validation probabilities."""
    precision, recall, thresholds = precision_recall_curve(y_true, proba)
    if len(thresholds) == 0:
        return 0.5, 0.0

    f1_scores = 2 * precision[1:] * recall[1:] / np.clip(precision[1:] + recall[1:], 1e-12, None)
    best_idx = int(np.nanargmax(f1_scores))
    return float(thresholds[best_idx]), float(f1_scores[best_idx])


def save_eval_plots(
    y_val: np.ndarray,
    val_proba: np.ndarray,
    y_test: np.ndarray,
    test_proba: np.ndarray,
    test_preds: np.ndarray,
    model_name: str,
    output_dir: str,
) -> List[str]:
    """Save ROC/PR/confusion plots and return artifact file paths."""
    plots_dir = Path(output_dir) / "figures"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # ROC + PR curves in one figure for quick model diagnostics.
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    fpr_v, tpr_v, _ = roc_curve(y_val, val_proba)
    fpr_t, tpr_t, _ = roc_curve(y_test, test_proba)
    auc_v = roc_auc_score(y_val, val_proba)
    auc_t = roc_auc_score(y_test, test_proba)

    axes[0].plot(fpr_v, tpr_v, label=f"Val AUC={auc_v:.3f}", linewidth=2)
    axes[0].plot(fpr_t, tpr_t, label=f"Test AUC={auc_t:.3f}", linewidth=2)
    axes[0].plot([0, 1], [0, 1], "k--", alpha=0.5)
    axes[0].set_title("ROC Curves")
    axes[0].set_xlabel("False Positive Rate")
    axes[0].set_ylabel("True Positive Rate")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    p_v, r_v, _ = precision_recall_curve(y_val, val_proba)
    p_t, r_t, _ = precision_recall_curve(y_test, test_proba)
    ap_v = average_precision_score(y_val, val_proba)
    ap_t = average_precision_score(y_test, test_proba)

    axes[1].plot(r_v, p_v, label=f"Val AP={ap_v:.3f}", linewidth=2)
    axes[1].plot(r_t, p_t, label=f"Test AP={ap_t:.3f}", linewidth=2)
    axes[1].set_title("Precision-Recall Curves")
    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    curve_path = plots_dir / f"{model_name.lower()}_roc_pr.png"
    fig.savefig(curve_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    cm = confusion_matrix(y_test, test_preds)
    fig_cm, ax_cm = plt.subplots(figsize=(5, 4))
    im = ax_cm.imshow(cm, cmap="Blues")
    ax_cm.set_title(f"{model_name} Test Confusion Matrix")
    ax_cm.set_xlabel("Predicted")
    ax_cm.set_ylabel("True")
    ax_cm.set_xticks([0, 1])
    ax_cm.set_yticks([0, 1])
    ax_cm.set_xticklabels(["non-binding", "binding"], rotation=15)
    ax_cm.set_yticklabels(["non-binding", "binding"])
    for (i, j), v in np.ndenumerate(cm):
        ax_cm.text(j, i, str(v), ha="center", va="center", color="black", fontsize=11)
    fig_cm.colorbar(im, ax=ax_cm, fraction=0.046, pad=0.04)
    plt.tight_layout()
    cm_path = plots_dir / f"{model_name.lower()}_test_confusion.png"
    fig_cm.savefig(cm_path, dpi=150, bbox_inches="tight")
    plt.close(fig_cm)

    return [str(curve_path), str(cm_path)]

def evaluate_model(
    model,
    X      : np.ndarray,
    y      : np.ndarray,
    threshold: float = 0.5,
    split  : str = "test",
) -> Tuple[dict, np.ndarray, np.ndarray]:
    """Return dict of metrics; also prints a classification report."""
    proba  = model.predict_proba(X)[:, 1]
    preds  = (proba >= threshold).astype(int)

    metrics = {
        f"{split}_auc"      : roc_auc_score(y, proba),
        f"{split}_pr_auc"   : average_precision_score(y, proba),
        f"{split}_accuracy" : accuracy_score(y, preds),
        f"{split}_precision": precision_score(y, preds, zero_division=0),
        f"{split}_recall"   : recall_score(y, preds, zero_division=0),
        f"{split}_threshold": float(threshold),
    }
    log.info(f"\n{'─'*50}\n{split.upper()} RESULTS\n{'─'*50}")
    for k, v in metrics.items():
        log.info(f"  {k:<25}: {v:.4f}")
    log.info("\n" + classification_report(y, preds, target_names=["non-binding", "binding"]))
    return metrics, proba, preds


# ─────────────────────────────────────────────────────────────────────────────
# 6. TRAINING WITH MLFLOW
# ─────────────────────────────────────────────────────────────────────────────

def train_model(
    name       : str,
    model,
    params     : dict,
    X_train    : np.ndarray,
    y_train    : np.ndarray,
    X_val      : np.ndarray,
    y_val      : np.ndarray,
    X_test     : np.ndarray,
    y_test     : np.ndarray,
    feat_names : List[str],
    cfg        : dict,
):
    """Train one model, evaluate, and log everything to MLflow."""
    with mlflow.start_run(run_name=name):
        # log config + model params
        mlflow.log_params({**params, "seed": cfg["seed"], "neg_ratio": cfg["neg_ratio"],
                           "use_kmers": cfg["use_kmers"], "no_esm": cfg["no_esm"]})

        log.info(f"\n{'═'*55}\nTraining: {name}\n{'═'*55}")
        model.fit(X_train, y_train)

        # Tune threshold on validation set using F1 to avoid fixed 0.5 cutoff.
        val_proba_for_thr = model.predict_proba(X_val)[:, 1]
        best_thr, best_val_f1 = choose_threshold_by_f1(y_val, val_proba_for_thr)
        mlflow.log_metric("val_best_f1", best_val_f1)
        mlflow.log_metric("val_best_threshold", best_thr)
        log.info(f"Chosen threshold from val (max F1): {best_thr:.4f} | F1={best_val_f1:.4f}")

        # evaluate on val + test
        val_metrics, val_proba, _ = evaluate_model(
            model, X_val, y_val, threshold=best_thr, split="val"
        )
        test_metrics, test_proba, test_preds = evaluate_model(
            model, X_test, y_test, threshold=best_thr, split="test"
        )

        # log all metrics
        mlflow.log_metrics({**val_metrics, **test_metrics})

        # log model artifact
        mlflow.sklearn.log_model(model, artifact_path=name)

        # log diagnostic plots as artifacts
        artifact_paths = save_eval_plots(
            y_val=y_val,
            val_proba=val_proba,
            y_test=y_test,
            test_proba=test_proba,
            test_preds=test_preds,
            model_name=name,
            output_dir=cfg["output_dir"],
        )
        for ap in artifact_paths:
            mlflow.log_artifact(ap, artifact_path="figures")

        # feature importances (RF only)
        clf = model.named_steps["clf"] if hasattr(model, "named_steps") else model
        if hasattr(clf, "feature_importances_"):
            imps = clf.feature_importances_
            top5 = np.argsort(imps)[-5:][::-1]
            log.info("Top 5 features (RF):")
            for i in top5:
                log.info(f"  {feat_names[i]:<35} {imps[i]:.4f}")

        return test_metrics


# ─────────────────────────────────────────────────────────────────────────────
# 7. MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def run_pipeline(cfg: dict):
    Path(cfg["output_dir"]).mkdir(parents=True, exist_ok=True)

    # ── load data ─────────────────────────────────────────────────────────────
    data_path = cfg["data_path"]
    if Path(data_path).exists():
        df_pos = load_from_anndata(data_path, cfg["max_samples"])
    else:
        raise FileNotFoundError(f"Data not found at {data_path}")

    # ── ESM2 embedder ─────────────────────────────────────────────────────────
    embedder = None
    if not cfg["no_esm"]:
        embedder = ESMEmbedder(cfg["esm_model"], batch_size=cfg["esm_batch"])

    # ── MLflow setup ──────────────────────────────────────────────────────────
    mlflow.set_tracking_uri(cfg["mlflow_uri"])
    mlflow.set_experiment(cfg["experiment"])

    results_by_model = {"LogisticRegression": [], "RandomForest": []}

    for i in range(cfg["n_seeds"]):
        run_seed = cfg["seed"] + i
        set_seed(run_seed)
        log.info(f"\n{'#' * 60}\nSEED RUN {i + 1}/{cfg['n_seeds']} | seed={run_seed}\n{'#' * 60}")

        # ── generate negatives ────────────────────────────────────────────────
        df = generate_negatives(df_pos, cfg["neg_ratio"], run_seed)

        # ── feature extraction ────────────────────────────────────────────────
        X, feat_names, _ = extract_features_batched(
            df, embedder, use_kmers=cfg["use_kmers"]
        )
        y = df["label"].values

        # ── split ─────────────────────────────────────────────────────────────
        train_idx, val_idx, test_idx = split_by_epitope(
            df, X, cfg["test_frac"], cfg["val_frac"], run_seed
        )
        X_train, y_train = X[train_idx], y[train_idx]
        X_val, y_val = X[val_idx], y[val_idx]
        X_test, y_test = X[test_idx], y[test_idx]

        # free raw feature matrix
        del X
        gc.collect()

        run_cfg = {**cfg, "seed": run_seed}

        # ── Logistic Regression ───────────────────────────────────────────────
        lr_params = {
            "model": "LogisticRegression",
            "C": cfg["lr_C"],
            "max_iter": cfg["lr_max_iter"],
        }
        lr_pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(
                C=cfg["lr_C"],
                max_iter=cfg["lr_max_iter"],
                class_weight="balanced",
                random_state=run_seed,
                solver="saga",
                n_jobs=-1,
            )),
        ])
        lr_metrics = train_model(
            "LogisticRegression", lr_pipe, lr_params,
            X_train, y_train, X_val, y_val, X_test, y_test, feat_names, run_cfg,
        )
        results_by_model["LogisticRegression"].append(lr_metrics)

        # ── Random Forest ─────────────────────────────────────────────────────
        rf_params = {
            "model": "RandomForest",
            "n_estimators": cfg["rf_n_trees"],
            "max_depth": cfg["rf_max_depth"],
        }
        rf_pipe = Pipeline([
            ("clf", RandomForestClassifier(
                n_estimators=cfg["rf_n_trees"],
                max_depth=cfg["rf_max_depth"],
                class_weight="balanced",
                random_state=run_seed,
                n_jobs=cfg["rf_n_jobs"],
            )),
        ])
        rf_metrics = train_model(
            "RandomForest", rf_pipe, rf_params,
            X_train, y_train, X_val, y_val, X_test, y_test, feat_names, run_cfg,
        )
        results_by_model["RandomForest"].append(rf_metrics)

    # ── Summary ───────────────────────────────────────────────────────────────
    log.info("\n" + "═" * 55)
    log.info("FINAL COMPARISON (test AUC / PR-AUC)")
    log.info("═" * 55)
    for model_name, metric_list in results_by_model.items():
        aucs = np.array([m["test_auc"] for m in metric_list], dtype=float)
        pr_aucs = np.array([m["test_pr_auc"] for m in metric_list], dtype=float)
        if len(metric_list) == 1:
            log.info(
                f"  {model_name:<25}: AUC = {aucs[0]:.4f} | PR-AUC = {pr_aucs[0]:.4f}"
            )
        else:
            log.info(
                f"  {model_name:<25}: "
                f"AUC = {aucs.mean():.4f} ± {aucs.std(ddof=1):.4f} | "
                f"PR-AUC = {pr_aucs.mean():.4f} ± {pr_aucs.std(ddof=1):.4f}"
            )
    log.info(f"\nMLflow UI: mlflow ui --backend-store-uri {cfg['mlflow_uri']}")


# ─────────────────────────────────────────────────────────────────────────────
# SMOKE TEST  — one batch only
# ─────────────────────────────────────────────────────────────────────────────

def one_batch_test(cfg: dict):
    """
    Runs the full pipeline on a tiny slice (first 200 positives) to verify
    every step works before submitting a long cluster job.
    """
    log.info("=== ONE-BATCH SMOKE TEST ===")
    cfg = {**cfg, "max_samples": 300, "neg_ratio": 1.0, "no_esm": True,
           "use_kmers": False}
    run_pipeline(cfg)
    log.info("Smoke test PASSED ✓")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="TCR–Epitope classical ML baseline")
    for k, v in DEFAULTS.items():
        if isinstance(v, bool):
            # support both --flag and --no_flag
            p.add_argument(f"--{k}", action="store_true",  default=v)
        elif v is None:
            p.add_argument(f"--{k}", type=str, default=v)
        else:
            p.add_argument(f"--{k}", type=type(v), default=v)
    return vars(p.parse_args())


if __name__ == "__main__":
    cfg = parse_args()
    log.info(f"Config: {cfg}")

    if cfg["mode"] == "one_batch":
        one_batch_test(cfg)
    elif cfg["mode"] == "train":
        run_pipeline(cfg)
    else:
        raise ValueError(f"Unknown mode '{cfg['mode']}'. Use 'one_batch' or 'train'.")