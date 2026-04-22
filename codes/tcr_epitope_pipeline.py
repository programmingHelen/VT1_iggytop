

"""
TCR–Epitope Binding — Classical ML Baseline
============================================
Mejoras respecto al código anterior:
  ✓ Generación de negativos (sin esto el modelo no aprende nada)
  ✓ StrictTCR split (sin data leakage)
  ✓ Cross-validation con AUC + AUPR
  ✓ Fácil de escalar cambiando los CONFIG params de arriba
 
Roadmap:
  Este archivo → ML baseline
  Siguiente     → GNN con PyTorch Geometric
"""
 
import os
import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmModel
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import StandardScaler
import anndata as ad
import warnings
warnings.filterwarnings("ignore")
 
 
# ─────────────────────────────────────────────
# CONFIG — cambia aquí para subir dificultad
# ─────────────────────────────────────────────
 
CONFIG = {
    # Cuántos pares positivos usar (None = todos)
    # Empieza con 200, luego sube a 1000, luego None
    "n_samples": 200,
 
    # Negativos por positivo. El paper GTE usa 10x por defecto
    # Empieza con 1, luego 5, luego 10
    "neg_ratio": 1,
 
    # Modelo ESM. Opciones por velocidad:
    #   "facebook/esm2_t6_8M_UR50D"   ← más rápido (8M params)
    #   "facebook/esm2_t12_35M_UR50D" ← mejor calidad
    #   "facebook/esm2_t33_650M_UR50D"← el que usa GTE (lento)
    "esm_model": "facebook/esm2_t6_8M_UR50D",
 
    # Batch size para embeddings
    "batch_size": 32,
 
    # Folds para cross-validation
    "n_folds": 5,
 
    # Random Forest params
    "rf_n_estimators": 200,
    "rf_max_depth": 10,
}
 
 
# ─────────────────────────────────────────────
# 1. LOAD DATA
# ─────────────────────────────────────────────
 
print("=" * 55)
print("TCR–Epitope ML Baseline")
print("=" * 55)
 
adata = ad.read_h5ad("data/deduplicated_anndata.h5ad")
 
 
# ─────────────────────────────────────────────
# 2. EXTRACT (cdr3b, epitope) — tu función, sin cambios
# ─────────────────────────────────────────────
 
def explode_tcr_epitope(adata):
    airr = adata.obsm["airr"]
    epitopes = adata.obs["epitope_sequence"]
    rows = []
    for i, chains in enumerate(airr):
        epitope = epitopes.iloc[i]
        for chain in chains:
            if chain is None:
                continue
            try:
                if chain["locus"] == "TRB" and chain["productive"]:
                    cdr3b = chain["junction_aa"]
                    if cdr3b is not None:
                        rows.append({"cdr3b": cdr3b, "epitope": epitope})
            except Exception:
                continue
    return pd.DataFrame(rows)
 
 
df_pos = explode_tcr_epitope(adata)
df_pos = df_pos.dropna()
df_pos = df_pos[df_pos["cdr3b"].str.len().between(8, 25)]
df_pos = df_pos.drop_duplicates(subset=["cdr3b", "epitope"])
df_pos["label"] = 1
 
print(f"\nPares positivos únicos: {len(df_pos)}")
print(f"Epítopos únicos: {df_pos['epitope'].nunique()}")
print(f"TCRs únicos: {df_pos['cdr3b'].nunique()}")
 
# Subsample si CONFIG["n_samples"] está definido
if CONFIG["n_samples"] and CONFIG["n_samples"] < len(df_pos):
    df_pos = df_pos.sample(CONFIG["n_samples"], random_state=42).reset_index(drop=True)
    print(f"→ Usando {len(df_pos)} pares (subsample)")
 
 
# ─────────────────────────────────────────────
# 3. GENERAR NEGATIVOS
#
# La clave: mezclar TCRs con epítopos que NO son su par cognato.
# Sin esto el modelo solo ve clase 1 y no aprende a distinguir.
# ─────────────────────────────────────────────
 
def generate_negatives(df_pos, neg_ratio=1, random_state=42):
    rng = np.random.default_rng(random_state)
    positive_pairs = set(zip(df_pos["cdr3b"], df_pos["epitope"]))
 
    tcrs = df_pos["cdr3b"].values
    epitopes = df_pos["epitope"].values
    n_neg = len(df_pos) * neg_ratio
 
    negatives = []
    attempts = 0
    max_attempts = n_neg * 20
 
    while len(negatives) < n_neg and attempts < max_attempts:
        tcr = rng.choice(tcrs)
        epi = rng.choice(epitopes)
        if (tcr, epi) not in positive_pairs:
            negatives.append({"cdr3b": tcr, "epitope": epi, "label": 0})
        attempts += 1
 
    df_neg = pd.DataFrame(negatives)
    df_full = pd.concat([df_pos, df_neg], ignore_index=True)
 
    print(f"\nPositivos: {len(df_pos)} | Negativos: {len(df_neg)} "
          f"(ratio 1:{neg_ratio})")
    return df_full
 
 
df_full = generate_negatives(df_pos, neg_ratio=CONFIG["neg_ratio"])
 
 
# ─────────────────────────────────────────────
# 4. ESM EMBEDDINGS con cache
#
# ESM2 convierte cada secuencia en un vector denso que captura
# propiedades evolutivas y estructurales — mucho mejor que one-hot.
# Mean pooling sobre los tokens → vector de tamaño fijo.
# ─────────────────────────────────────────────
 
print(f"\nCargando ESM: {CONFIG['esm_model']}...")
tokenizer = AutoTokenizer.from_pretrained(CONFIG["esm_model"])
esm = EsmModel.from_pretrained(CONFIG["esm_model"])
esm.eval()
device = "cuda" if torch.cuda.is_available() else "cpu"
esm = esm.to(device)
print(f"Dispositivo: {device}")
 
embedding_cache = {}
 
def embed_sequences(seqs, batch_size=16):
    unique_seqs = list(dict.fromkeys(seqs))  # deduplica preservando orden
    to_embed = [s for s in unique_seqs if s not in embedding_cache]
 
    for i in range(0, len(to_embed), batch_size):
        batch = to_embed[i:i + batch_size]
        inputs = tokenizer(batch, return_tensors="pt", padding=True,
                           truncation=True, max_length=50).to(device)
        with torch.no_grad():
            out = esm(**inputs)
        emb = out.last_hidden_state.mean(dim=1).cpu().numpy()
        for s, e in zip(batch, emb):
            embedding_cache[s] = e
 
    return np.array([embedding_cache[s] for s in seqs])
 
 
print("Embedding secuencias únicas (TCR + Epitopo)...")
all_seqs = list(set(df_full["cdr3b"].tolist() + df_full["epitope"].tolist()))
embed_sequences(all_seqs, batch_size=CONFIG["batch_size"])
print(f"  Cache: {len(embedding_cache)} secuencias embedidas")
 
# Construir X: concatenar embedding TCR + embedding Epitopo
tcr_emb = np.array([embedding_cache[s] for s in df_full["cdr3b"]])
epi_emb = np.array([embedding_cache[s] for s in df_full["epitope"]])
X = np.concatenate([tcr_emb, epi_emb], axis=1)
y = df_full["label"].values
 
print(f"  X shape: {X.shape}  →  {X.shape[1]} features por par")
 
 
# ─────────────────────────────────────────────
# 5. STRICT TCR SPLIT
#
# Problema con split aleatorio: el mismo TCR puede aparecer
# en train y en test → el modelo "recuerda" el TCR, no generaliza.
# StrictTCR garantiza que los TCRs del test nunca se vieron en train.
# Esto es más realista (predicción sobre TCRs nuevos).
# ─────────────────────────────────────────────
 
def strict_tcr_kfold(df, n_splits=5, random_state=42):
    """
    Divide por TCR único, no por par.
    Garantiza que ningún TCR del test esté en train.
    """
    unique_tcrs = df["cdr3b"].unique()
    rng = np.random.default_rng(random_state)
    rng.shuffle(unique_tcrs)
 
    folds = []
    tcr_chunks = np.array_split(unique_tcrs, n_splits)
 
    for i in range(n_splits):
        test_tcrs = set(tcr_chunks[i])
        test_idx = df.index[df["cdr3b"].isin(test_tcrs)].tolist()
        train_idx = df.index[~df["cdr3b"].isin(test_tcrs)].tolist()
        folds.append((train_idx, test_idx))
 
    return folds
 
 
# ─────────────────────────────────────────────
# 6. CROSS-VALIDATION
# ─────────────────────────────────────────────
 
print(f"\nCross-validation ({CONFIG['n_folds']}-fold StrictTCR)...")
print("─" * 45)
 
folds = strict_tcr_kfold(df_full, n_splits=CONFIG["n_folds"])
scaler = StandardScaler()
 
auc_scores, aupr_scores = [], []
 
for fold, (train_idx, test_idx) in enumerate(folds):
    if len(test_idx) == 0 or len(set(y[test_idx])) < 2:
        print(f"  Fold {fold+1}: skip (sin ambas clases en test)")
        continue
 
    X_train = scaler.fit_transform(X[train_idx])
    X_test = scaler.transform(X[test_idx])
    y_train, y_test = y[train_idx], y[test_idx]
 
    clf = RandomForestClassifier(
        n_estimators=CONFIG["rf_n_estimators"],
        max_depth=CONFIG["rf_max_depth"],
        class_weight="balanced",   # compensa desbalance pos/neg
        n_jobs=-1,
        random_state=42,
    )
    clf.fit(X_train, y_train)
    probs = clf.predict_proba(X_test)[:, 1]
 
    auc = roc_auc_score(y_test, probs)
    aupr = average_precision_score(y_test, probs)
    auc_scores.append(auc)
    aupr_scores.append(aupr)
    print(f"  Fold {fold+1}: AUC={auc:.3f}  AUPR={aupr:.3f}  "
          f"(train={len(train_idx)}, test={len(test_idx)})")
 
print("─" * 45)
if auc_scores:
    print(f"  AUC  medio: {np.mean(auc_scores):.3f} ± {np.std(auc_scores):.3f}")
    print(f"  AUPR medio: {np.mean(aupr_scores):.3f} ± {np.std(aupr_scores):.3f}")
    print("─" * 45)
 
print("""
Referencia paper GTE (StrictTCR, ESM2-650M):
  AUC  ~0.82  AUPR ~0.47  (VDJdb)
  AUC  ~0.90  AUPR ~0.63  (pMTnet)
 
Próximos pasos para mejorar este baseline:
  1. Subir neg_ratio a 5 o 10
  2. Subir n_samples a None (todos los datos)
  3. Cambiar esm_model a esm2_t12_35M o esm2_t33_650M
  4. → GNN con PyTorch Geometric
""")