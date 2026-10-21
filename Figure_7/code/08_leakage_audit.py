"""
08_leakage_audit.py
---------------------
Answers the editor's concern with NUMBERS: how optimistic is a random
row/chromosome split compared with the organism-level split?

For every Group2 with enough organisms, the same RandomForest (100 trees,
same inputs as the main ML analysis) is evaluated twice with 5-fold CV:

  RF_random_rows     : rows (chromosomes) split at random -> sibling chromosomes
                       of every test chromosome are in the training set (leaky)
  RF_organism_level  : StratifiedGroupKFold by organism (the strategy used by the
                       main analysis)

Scores are Pearson r (all chromosomes), MAE, and Within_org_r (can the model predict
how a chromosome differs from its OWN organism's mean?). Pearson r falling and MAE
rising when going from the random split to the organism-level split = optimism that a
random split would have hidden; report the table in the Supplement. Within_org_r is a
separate diagnostic and does NOT necessarily move in the same direction (on the
synthetic test data it went UP under the organism-level split), so interpret it on
its own and do not use it as a leakage indicator.

Output: tables/leakage_audit.csv (+ figure LeakageAudit if matplotlib is available)
"""
import importlib
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold

warnings.filterwarnings("ignore")
ml_mod = importlib.import_module("03_ml_methods")
scoring = importlib.import_module("scoring")
gs = importlib.import_module("grouped_split")
BASES = ml_mod.BASES


def _rf_fold(task):
    g, scheme, fold, Xtr, Ytr, Xte, Yte, rid, uc = task
    out = {}
    for j, b in enumerate(BASES):
        m = RandomForestRegressor(n_estimators=100, random_state=42 + fold, n_jobs=1).fit(Xtr, Ytr[:, j])
        out[b] = m.predict(Xte)
    return g, scheme, fold, rid, uc, Yte, out


def run_leakage_audit(df_eval, valid_groups, n_splits=5, n_jobs=-1, min_organisms=10,
                      seed=42, verbose=True):
    t0 = time.time()
    tasks = []
    for g in valid_groups:
        dfg = df_eval[df_eval["Group2"] == g]
        if len(dfg) == 0:
            continue
        X, Y, meta = ml_mod._build_group_features(dfg)
        if len(X) < ml_mod.MIN_SAMPLES or meta["SplitGroup"].nunique() < min_organisms:
            continue
        Xa, Ya = X.to_numpy(), Y.to_numpy()
        rid, uc = meta["RowID"].to_numpy(), meta["UnitCode"].to_numpy()
        org_folds = [(f, tr, te) for _, f, tr, te in gs.grouped_folds(
            meta, group_col="SplitGroup", n_splits=n_splits, n_repeats=1, seed=seed, min_groups=min_organisms)]
        rnd_folds = [(f, tr, te) for f, (tr, te) in enumerate(
            KFold(n_splits, shuffle=True, random_state=seed).split(Xa))]
        for scheme, folds in (("RF_organism_level", org_folds), ("RF_random_rows", rnd_folds)):
            for f, tr, te in folds:
                tasks.append((g, scheme, f, Xa[tr], Ya[tr], Xa[te], Ya[te], rid[te], uc[te]))
    if not tasks:
        print("[WARN] leakage audit: no group has enough organisms.")
        return pd.DataFrame()
    tasks.sort(key=lambda t: -len(t[3]))
    if verbose:
        print(f"[INFO] Leakage audit: {len(tasks)} RandomForest fold-tasks ...", flush=True)
    res = Parallel(n_jobs=n_jobs)(delayed(_rf_fold)(t) for t in tasks)
    frames = []
    for g, scheme, f, rid, uc, Yte, out in res:
        for j, b in enumerate(BASES):
            frames.append(pd.DataFrame({"RowID": rid, "UnitCode": uc, "Group2": g, "Method": scheme,
                                        "Base": b, "Repeat": 0, "Fold": f,
                                        "Predicted": out[b], "Observed": Yte[:, j]}))
    long = pd.concat(frames, ignore_index=True)
    sc = scoring.score_per_group(long, n_boot=0)
    tab = sc.pivot_table(index=["Group2", "Base"], columns="Method",
                         values=["Pearson_r", "Within_org_r", "MAE"]).reset_index()
    tab.columns = ["_".join([c for c in col if c]) for col in tab.columns.to_flat_index()]
    for m in ("Pearson_r", "Within_org_r"):
        tab[f"{m}_drop"] = tab[f"{m}_RF_random_rows"] - tab[f"{m}_RF_organism_level"]
    tab["MAE_increase"] = tab["MAE_RF_organism_level"] - tab["MAE_RF_random_rows"]
    if verbose:
        print(f"[INFO] Leakage audit done in {time.time() - t0:.1f}s")
    return tab


def summarize(tab):
    d = tab.groupby("Group2")[["Pearson_r_RF_random_rows", "Pearson_r_RF_organism_level", "Pearson_r_drop",
                               "Within_org_r_RF_random_rows", "Within_org_r_RF_organism_level",
                               "Within_org_r_drop", "MAE_increase"]].mean().round(3)
    return d


if __name__ == "__main__":
    paths = importlib.import_module("paths")
    df_eval = pd.read_csv(paths.table_path("df_eval.csv"))
    load_validate = importlib.import_module("01_load_validate")
    groups = [g for g in load_validate.GROUP_ORDER if g in set(df_eval["Group2"])]
    tab = run_leakage_audit(df_eval, groups)
    tab.to_csv(paths.table_path("leakage_audit.csv"), index=False)
    print(summarize(tab).to_string())
