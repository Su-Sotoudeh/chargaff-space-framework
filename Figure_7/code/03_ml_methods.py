"""
03_ml_methods.py
------------------
Trains 4 ML regressors (RandomForest, GradientBoost, SVR, NeuralNet) to
predict chromosome nucleotide composition (A_W, T_W, G_W, C_W) from the
annotation composition of BOTH strands + Z (FEATURE_SET = "plus_and_comp",
the same information that the dominant-strand formula M4 uses, so the
comparison isolates modelling approach rather than information advantage).

DESIGN (v2 -- organism-level evaluation)
----------------------------------------
1. PER-GROUP independence (unchanged): every Group2 gets its own models; a
   RandomForest fitted on Chordata never sees a Plasmid row.

2. ORGANISM-LEVEL SPLIT (new). Inside each Group2 the unit that is split is
   the ORGANISM (Organ_ID), never the chromosome or the row. All chromosomes
   of an organism are either all in training or all in testing. This removes
   the leakage criticised by the editor (train on chr1, test on chr2 of the
   same organism: they share GC content, Z, annotation composition ...).
     * StratifiedGroupKFold, 5 folds, strata = organism-mean GC quartile, so
       every test fold spans low- and high-GC organisms;
     * repeated N_REPEATS times with different seeds (default 3) -> mean and
       SD across repeats;
     * every row gets exactly ONE out-of-fold prediction per repeat, so the
       formula methods can be scored on exactly the same rows.
   Set SPLIT_GROUP_COL to a coarser taxonomy column (genus / family) to run a
   stricter "leave-clade-out" version.

3. Groups with fewer than MIN_ORGANISMS organisms (default 10) or fewer than
   MIN_SAMPLES rows are skipped for ML: with e.g. 6 organisms an organism-level
   5-fold split has no statistical meaning. Formula methods are still scored.

4. Folds are fitted in PARALLEL (joblib, all cores by default). Models inside
   a task are single-threaded to avoid oversubscription.

Notes
  * MLPRegressor(early_stopping=True) holds out a random 10% of the TRAINING
    fold for early stopping only; this never touches the test organisms.
  * SVR subsamples to SVR_MAX_TRAIN_ROWS_PER_GROUP training rows per fit.
  * No hyperparameter tuning is done (defaults below), so there is no
    selection on the test organisms. If you add tuning, nest it (inner
    GroupKFold on the training organisms only).

Output (long table): RowID, UnitCode, Group2, Method, Base, Repeat, Fold,
Predicted, Observed, Split("test" = prediction made by a model that never saw
this organism).
"""

import importlib
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

warnings.filterwarnings("ignore")
scoring = importlib.import_module("scoring")
gs = importlib.import_module("grouped_split")
ORGAN_COL = gs.ORGAN_COL

BASES = ["A", "T", "G", "C"]
ML_METHODS = ["RandomForest", "GradientBoost", "SVR", "NeuralNet"]
RANDOM_STATE = 42

# ---- evaluation design ----
SPLIT_GROUP_COL = ORGAN_COL   # or "Genus"/"Family" if you add such a column
N_SPLITS = 5
N_REPEATS = 3
MIN_ORGANISMS = 10            # fewer organisms than this -> group skipped for ML
MIN_SAMPLES = 30              # fewer rows than this -> group skipped for ML
N_JOBS = -1                   # parallel workers (-1 = all cores)

# ---- model inputs ----
# "plus_and_comp" (default) gives ML both annotated strands + Z, i.e. at least the
# information the dominant-strand formula M4 uses, and works for single-stranded
# rows (Z = 0 / 100). "plus_only" gives ML only the +annotation and Z.
FEATURE_SETS = {
    "plus_only": {"annot_A": "annot_P_A", "annot_T": "annot_P_T", "annot_G": "annot_P_G",
                  "annot_C": "annot_P_C", "Z": "Z"},
    # more information than Formula M1 (M2/M4/M5 use the complement annotation):
    "plus_and_comp": {"annot_A": "annot_P_A", "annot_T": "annot_P_T", "annot_G": "annot_P_G",
                      "annot_C": "annot_P_C", "Z": "Z", "cA": "annot_C_A", "cT": "annot_C_T",
                      "cG": "annot_C_G", "cC": "annot_C_C", "C_portion": "C_portion"},
}
FEATURE_SET = "plus_and_comp"   # same information as M4 (both strands); "plus_only" = M1 information only

SVR_MAX_TRAIN_ROWS_PER_GROUP = 5000


def make_model(name, seed):
    """One fresh single-threaded model (parallelism is across folds)."""
    if name == "RandomForest":
        return RandomForestRegressor(n_estimators=200, random_state=seed, n_jobs=1)
    if name == "GradientBoost":
        return HistGradientBoostingRegressor(max_iter=200, max_depth=3, random_state=seed)
    if name == "SVR":
        return SVR(kernel="rbf", C=10, epsilon=0.1)
    if name == "NeuralNet":
        return MLPRegressor(hidden_layer_sizes=(50, 50), activation="relu", solver="adam",
                            max_iter=1000, random_state=seed, early_stopping=True)
    raise ValueError(name)


def _build_group_features(dfg, feature_set=FEATURE_SET):
    """X, Y and a meta frame aligned row-for-row. df_eval is already
    restricted to finite/evaluable rows (see 02_formula_methods)."""
    cols = FEATURE_SETS[feature_set]
    X = pd.DataFrame({new: pd.to_numeric(dfg[old], errors="coerce").to_numpy(float)
                      for new, old in cols.items()})
    Y = pd.DataFrame({b: pd.to_numeric(dfg[f"{b}_W"], errors="coerce").to_numpy(float) for b in BASES})
    if not (np.isfinite(X.to_numpy()).all() and np.isfinite(Y.to_numpy()).all()):
        raise ValueError("Non-finite values reached the ML stage; run restrict_to_evaluable() first.")
    split_col = dfg[SPLIT_GROUP_COL] if SPLIT_GROUP_COL in dfg.columns else dfg[ORGAN_COL]
    split_col = split_col.where(split_col.notna(), dfg[ORGAN_COL]).astype(str).to_numpy()
    meta = pd.DataFrame({
        "RowID": dfg["RowID"].to_numpy(), "UnitCode": dfg["UnitCode"].to_numpy(),
        "SplitGroup": split_col, "G_W": Y["G"].to_numpy(), "C_W": Y["C"].to_numpy()})
    return X, Y, meta


def _fit_predict_fold(task):
    """Worker: fit all 4 models x 4 bases on one training fold, predict the test fold."""
    (g, rep, fold, Xtr, Ytr, Xte, Yte, rid_te, uc_te) = task
    seed = RANDOM_STATE + 1000 * rep + fold
    scaler = StandardScaler().fit(Xtr)
    Xtr_s, Xte_s = scaler.transform(Xtr), scaler.transform(Xte)
    if len(Xtr) > SVR_MAX_TRAIN_ROWS_PER_GROUP:
        svr_idx = np.random.RandomState(seed).choice(len(Xtr), SVR_MAX_TRAIN_ROWS_PER_GROUP, replace=False)
    else:
        svr_idx = slice(None)
    preds = {}
    for name in ML_METHODS:
        for j, base in enumerate(BASES):
            model = make_model(name, seed)
            if name == "SVR":
                model.fit(Xtr_s[svr_idx], Ytr[svr_idx, j]); p = model.predict(Xte_s)
            elif name == "NeuralNet":
                model.fit(Xtr_s, Ytr[:, j]); p = model.predict(Xte_s)
            else:
                model.fit(Xtr, Ytr[:, j]); p = model.predict(Xte)
            preds[(name, base)] = p
    return g, rep, fold, rid_te, uc_te, Yte, preds


def run_ml_methods_per_group(df_eval, valid_groups, n_splits=N_SPLITS, n_repeats=N_REPEATS,
                             min_organisms=MIN_ORGANISMS, min_samples=MIN_SAMPLES,
                             n_jobs=N_JOBS, feature_set=FEATURE_SET, verbose=True):
    """Organism-level repeated CV, independently for every Group2.
    Returns (ml_results_long, skipped_groups[(group, n_rows, n_organisms)])."""
    t0 = time.time()
    tasks, skipped, order = [], [], []
    for g in valid_groups:
        dfg = df_eval[df_eval["Group2"] == g]
        if len(dfg) == 0:
            continue
        X, Y, meta = _build_group_features(dfg, feature_set)
        n_org = meta["SplitGroup"].nunique()
        if len(X) < min_samples or n_org < min_organisms:
            skipped.append((g, len(X), n_org))
            continue
        folds = gs.grouped_folds(meta, group_col="SplitGroup", n_splits=n_splits,
                                 n_repeats=n_repeats, seed=RANDOM_STATE, min_groups=min_organisms)
        if verbose:
            print(f"[INFO] Group '{g}': {len(X)} rows, {n_org} organisms -> "
                  f"{len(folds)} organism-level folds ({n_repeats} repeats x {min(n_splits, n_org)})", flush=True)
        Xa, Ya = X.to_numpy(), Y.to_numpy()
        for rep, fold, tr, te in folds:
            tasks.append((g, rep, fold, Xa[tr], Ya[tr], Xa[te], Ya[te],
                          meta["RowID"].to_numpy()[te], meta["UnitCode"].to_numpy()[te]))
    # biggest tasks first -> better load balance
    tasks.sort(key=lambda t: -len(t[3]))
    if skipped and verbose:
        print(f"[WARN] {len(skipped)} group(s) skipped for ML (< {min_organisms} organisms or "
              f"< {min_samples} rows) -- (group, rows, organisms): {skipped}")
    if not tasks:
        return pd.DataFrame(), skipped
    if verbose:
        print(f"[INFO] Running {len(tasks)} fold-tasks in parallel (n_jobs={n_jobs}) ...", flush=True)
    results = Parallel(n_jobs=n_jobs, verbose=5 if verbose else 0)(
        delayed(_fit_predict_fold)(t) for t in tasks)

    frames = []
    for g, rep, fold, rid, uc, Yte, preds in results:
        for (name, base), p in preds.items():
            j = BASES.index(base)
            frames.append(pd.DataFrame({
                "RowID": rid, "UnitCode": uc, "Group2": g, "Method": name, "Base": base,
                "Repeat": np.int8(rep), "Fold": np.int8(fold), "Predicted": p, "Observed": Yte[:, j],
                "Split": "test"}))
    ml_long = pd.concat(frames, ignore_index=True)
    for c in ("Group2", "Method", "Base", "Split"):
        ml_long[c] = ml_long[c].astype("category")
    if verbose:
        print(f"[INFO] ML done in {time.time() - t0:.1f}s", flush=True)
    return ml_long, skipped


def score_ml_methods_per_group(ml_results_long, min_samples=3, n_boot=1000):
    scores = scoring.score_per_group(ml_results_long, n_boot=n_boot, min_samples=min_samples)
    return scores, scoring.best_per_group_dict(scores)


def score_ml_methods_pooled(ml_results_long, n_boot=200):
    scores = scoring.score_pooled(ml_results_long, n_boot=n_boot)
    return scores, scoring.best_pooled_dict(scores)


if __name__ == "__main__":
    paths = importlib.import_module("paths")
    load_validate = importlib.import_module("01_load_validate")
    formula_mod = importlib.import_module("02_formula_methods")

    df_clean = pd.read_csv(paths.table_path("df_clean.csv"))
    valid_groups = [g for g in load_validate.GROUP_ORDER if g in set(df_clean["Group2"])]
    df_eval, _, valid_eval = formula_mod.restrict_to_evaluable(df_clean, valid_groups, load_validate.GROUP_ORDER)

    ml_long, skipped = run_ml_methods_per_group(df_eval, valid_eval)
    scores_per_group, _ = score_ml_methods_per_group(ml_long)
    scores_pooled, best_pooled = score_ml_methods_pooled(ml_long)
    ml_long.to_csv(paths.table_path("ml_results_long.csv"), index=False)
    scores_per_group.to_csv(paths.table_path("ml_scores_per_group.csv"), index=False)
    scores_pooled.to_csv(paths.table_path("ml_scores_pooled.csv"), index=False)
    print(scores_pooled.sort_values(["Base", "Pearson_r"], ascending=[True, False])
          [["Method", "Base", "Pearson_r", "Pearson_r_lo", "Pearson_r_hi", "CCC", "MAE", "N"]].to_string(index=False))
