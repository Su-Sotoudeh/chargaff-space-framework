"""
grouped_split.py
-----------------
Organism-level (leakage-free) cross-validation helpers.

WHY: chromosomes of one organism share GC content, codon bias, Z, etc.
Splitting rows (or chromosomes) at random puts siblings of every test
chromosome in the training set, so the test is optimistic. Here an organism
(Organ_ID) is NEVER split across train and test.

DESIGN
  * StratifiedGroupKFold: groups = Organ_ID, strata = organism-mean GC
    quartile -> every test fold spans low/high GC organisms.
  * Every row receives exactly one out-of-fold (OOF) prediction per repeat,
    so formula methods (no training) are scored on EXACTLY the same rows.
  * Repeats with different seeds -> mean +/- SD across repeats.
  * SPLIT_GROUP_COL can be set to a coarser taxonomy column (e.g. "Genus",
    "Family") if you add one to the CSV -> "leave-clade-out" robustness test.
"""
import warnings
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

# Column names in your CSV (change here only)
ORGAN_COL = "Organ_ID"   # organism id
CHR_COL = "Chr_ID"       # chromosome id within organism


def organism_gc_strata(meta, group_col=ORGAN_COL, n_bins=4):
    """One stratum label per organism (GC quantile), broadcast to its rows."""
    gc = (meta["G_W"] + meta["C_W"]).groupby(meta[group_col]).mean()
    n_bins = int(min(n_bins, max(2, len(gc))))
    bins = pd.qcut(gc.rank(method="first"), q=n_bins, labels=False)
    return meta[group_col].map(bins).to_numpy()


def grouped_folds(meta, group_col=ORGAN_COL, n_splits=5, n_repeats=3,
                  seed=42, min_groups=10, stratify=True):
    """
    meta: DataFrame aligned row-for-row with X/Y, containing group_col and
          (if stratify) G_W, C_W.
    Returns list of (repeat, fold, train_pos, test_pos) with POSITIONAL indices.
    Returns [] if there are fewer than `min_groups` distinct groups (skip it).
    """
    groups = meta[group_col].to_numpy()
    n_org = len(np.unique(groups))
    if n_org < min_groups:
        return []
    k = min(n_splits, n_org)
    y = organism_gc_strata(meta, group_col) if stratify else np.zeros(len(meta), int)
    dummy_X = np.zeros(len(meta))
    out = []
    for rep in range(n_repeats):
        cv = StratifiedGroupKFold(n_splits=k, shuffle=True, random_state=seed + rep)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for fold, (tr, te) in enumerate(cv.split(dummy_X, y, groups)):
                if not set(groups[tr]).isdisjoint(groups[te]):
                    raise AssertionError("group leaked across train/test split!")
                out.append((rep, fold, tr, te))
    return out
