"""
02_formula_methods.py
----------------------
The 5 formula-based prediction methods (M1-M5), plus a baseline (M0).

  M0_Identity : the composition of the DOMINANT annotated strand (the same
                strand M4 uses), normalised to 100% and used as-is as the
                whole-chromosome composition -- i.e. M4 WITHOUT the Chargaff
                correction. This is the "do-nothing" baseline: if M1-M5 and
                the ML models do not beat it clearly (see MAE_skill_vs_M0), a
                high Pearson r only reflects that annotation GC tracks
                chromosome GC.

COMPARABILITY: the headline comparison (best formula M4 vs ML vs M0) is made
on exactly the same rows. `restrict_to_evaluable()` keeps every row for which
the observed composition, M0 and M4 are finite. Single-stranded groups
(Z = 0 or Z = 100: ssDNA/ssRNA viruses, many plasmids) are therefore KEPT:
M4 and M0 simply use the one strand that carries annotation.
M1, M2, M3 and M5 each need one specific strand and are undefined where that
strand is empty (e.g. M1 at Z = 0, M2 at Z = 100). They are scored on the rows
where they are defined (the N column shows how many); their skill vs M0 is
computed on those same rows.

Formula methods use no training, so there is no train/test leakage here; the
organism-level split only matters for ML (03_ml_methods.py). Scores for the
formula methods are computed with the same code (scoring.py) as the ML scores,
including organism-bootstrap confidence intervals.

Output tables (long format): RowID, UnitCode, Group2, Method, Base, Repeat,
Fold, Predicted, Observed, Split.
"""

import importlib
import numpy as np
import pandas as pd

scoring = importlib.import_module("scoring")
gs = importlib.import_module("grouped_split")
ORGAN_COL, CHR_COL = gs.ORGAN_COL, gs.CHR_COL

BASES = ["A", "T", "G", "C"]
FORMULA_METHODS = ["M1_PlusOnly", "M2_CompOnly", "M3_Average", "M4_Dominant", "M5_Weighted"]
BASELINE_METHOD = "M0_Identity"


def compute_formula_prediction(A_plus, T_plus, G_plus, C_plus, Z):
    """
    Core Chargaff-style formula:
    given plus-strand (or comp-strand) annotated composition and Z
    (the percentage weight of that strand's annotation), reconstruct
    the predicted minus-strand contribution and combine.
    """
    A_plus = np.asarray(A_plus, dtype=float)
    T_plus = np.asarray(T_plus, dtype=float)
    G_plus = np.asarray(G_plus, dtype=float)
    C_plus = np.asarray(C_plus, dtype=float)
    Z = np.asarray(Z, dtype=float)

    with np.errstate(divide='ignore', invalid='ignore'):
        A_minus = T_plus * (100 - Z) / Z
        T_minus = A_plus * (100 - Z) / Z
        G_minus = C_plus * (100 - Z) / Z
        C_minus = G_plus * (100 - Z) / Z

        A = A_plus + A_minus
        T = T_plus + T_minus
        G = G_plus + G_minus
        C = C_plus + C_minus

        total = A + T + G + C

        A = 100 * A / total
        T = 100 * T / total
        G = 100 * G / total
        C = 100 * C / total

    return A, T, G, C


def _num(df, col):
    return pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)


def predict_all_methods(dfg):
    """Vectorised predictions of M0..M5 for the rows of dfg.
    Returns {method: (A, T, G, C) numpy arrays}. Non-finite = not evaluable."""
    Z = _num(dfg, "Z")
    P_portion = _num(dfg, "P_portion")
    C_portion = _num(dfg, "C_portion")
    Ap, Tp, Gp, Cp = (_num(dfg, f"annot_P_{b}") for b in BASES)
    Ac, Tc, Gc, Cc = (_num(dfg, f"annot_C_{b}") for b in BASES)

    out = {}
    out["M1_PlusOnly"] = compute_formula_prediction(Ap, Tp, Gp, Cp, Z)
    out["M2_CompOnly"] = compute_formula_prediction(Ac, Tc, Gc, Cc, C_portion)

    with np.errstate(divide="ignore", invalid="ignore"):
        total = Ap + Tp + Gp + Cp
        An, Tn, Gn, Cn = (100 * v / total for v in (Ap, Tp, Gp, Cp))
    out["M3_Average"] = ((An + Tn) / 2, (An + Tn) / 2, (Gn + Cn) / 2, (Gn + Cn) / 2)

    use_P = P_portion > C_portion
    A_dom, T_dom = np.where(use_P, Ap, Ac), np.where(use_P, Tp, Tc)
    G_dom, C_dom = np.where(use_P, Gp, Gc), np.where(use_P, Cp, Cc)
    Z_dom = np.where(use_P, P_portion, C_portion)
    out["M4_Dominant"] = compute_formula_prediction(A_dom, T_dom, G_dom, C_dom, Z_dom)

    A5 = (Ap * Z + Ac * (100 - Z)) / 100
    T5 = (Tp * Z + Tc * (100 - Z)) / 100
    G5 = (Gp * Z + Gc * (100 - Z)) / 100
    C5 = (Cp * Z + Cc * (100 - Z)) / 100
    out["M5_Weighted"] = compute_formula_prediction(A5, T5, G5, C5, Z)

    with np.errstate(divide="ignore", invalid="ignore"):
        tot_dom = A_dom + T_dom + G_dom + C_dom
        out[BASELINE_METHOD] = tuple(100 * v / tot_dom for v in (A_dom, T_dom, G_dom, C_dom))
    return out


CORE_METHODS = [BASELINE_METHOD, "M4_Dominant"]   # must be finite for a row to be evaluated


def evaluable_mask(df):
    """True for rows where the observed composition AND the core methods
    (M0 baseline, M4 dominant-strand formula) are finite for all 4 bases.
    These rows are used for ML as well, so the headline comparison is on
    identical rows. M1/M2/M3/M5 may still be undefined on some of them."""
    preds = predict_all_methods(df)
    ok = np.ones(len(df), dtype=bool)
    for m in CORE_METHODS:
        for a in preds[m]:
            ok &= np.isfinite(a)
    for c in ("A_W", "T_W", "G_W", "C_W"):
        ok &= np.isfinite(_num(df, c))
    return ok


def restrict_to_evaluable(df_clean, valid_groups, group_order=None):
    """
    Returns (df_eval, exclusion_table, valid_groups_eval).
    df_eval has RowID (0..n-1, unique) and UnitCode (integer id of
    (Group2, Organ_ID) = the independent unit used for splitting/bootstrapping).
    """
    df = df_clean.reset_index(drop=True).copy()
    ok = evaluable_mask(df)
    rows = []
    for g in valid_groups:
        m = (df["Group2"] == g).to_numpy()
        rows.append({"Group2": g, "n_rows": int(m.sum()), "n_excluded": int((m & ~ok).sum()),
                     "pct_excluded": round(100 * (m & ~ok).sum() / max(m.sum(), 1), 2)})
    excl = pd.DataFrame(rows)
    df_eval = df[ok].reset_index(drop=True)
    df_eval["RowID"] = np.arange(len(df_eval), dtype=np.int64)
    key = df_eval["Group2"].astype(str) + "||" + df_eval[ORGAN_COL].astype(str)
    df_eval["UnitCode"] = pd.factorize(key)[0].astype(np.int64)
    counts = df_eval["Group2"].value_counts()
    order = group_order or valid_groups
    valid_eval = [g for g in order if g in valid_groups and counts.get(g, 0) >= 3]
    return df_eval, excl, valid_eval


def _long_for(df_eval, valid_groups, methods, min_samples=3):
    frames = []
    for g in valid_groups:
        dfg = df_eval[df_eval["Group2"] == g]
        if len(dfg) < min_samples:
            continue
        preds = predict_all_methods(dfg)
        obs = {b: _num(dfg, f"{b}_W") for b in BASES}
        rid, uc = dfg["RowID"].to_numpy(), dfg["UnitCode"].to_numpy()
        for m in methods:
            # a method is scored only where it is defined for ALL 4 bases
            fin = np.ones(len(dfg), dtype=bool)
            for a in preds[m]:
                fin &= np.isfinite(a)
            if fin.sum() < min_samples:
                continue
            for i, b in enumerate(BASES):
                frames.append(pd.DataFrame({
                    "RowID": rid[fin], "UnitCode": uc[fin], "Group2": g, "Method": m, "Base": b,
                    "Predicted": preds[m][i][fin], "Observed": obs[b][fin]}))
    long = pd.concat(frames, ignore_index=True)
    long["Repeat"] = np.int8(0)
    long["Fold"] = np.int8(-1)
    long["Split"] = "all"
    for c in ("Group2", "Method", "Base", "Split"):
        long[c] = long[c].astype("category")
    return long[["RowID", "UnitCode", "Group2", "Method", "Base", "Repeat", "Fold",
                 "Predicted", "Observed", "Split"]]


def run_formula_methods(df_eval, valid_groups, min_samples=3):
    """M1-M5 for every group INDEPENDENTLY, on df_eval (already evaluable)."""
    return _long_for(df_eval, valid_groups, FORMULA_METHODS, min_samples)


def run_baseline_identity(df_eval, valid_groups, min_samples=3):
    """M0 baseline, same rows as everything else."""
    return _long_for(df_eval, valid_groups, [BASELINE_METHOD], min_samples)


def score_formula_methods_per_group(results_long, min_samples=3, n_boot=1000):
    """PRIMARY table: per (Group2, Method, Base). Returns (scores, best_dict)."""
    scores = scoring.score_per_group(results_long, n_boot=n_boot, min_samples=min_samples)
    return scores, scoring.best_per_group_dict(scores)


def score_formula_methods_pooled(results_long, n_boot=200):
    """SECONDARY headline table (dominated by big groups). Returns (scores, best_dict)."""
    scores = scoring.score_pooled(results_long, n_boot=n_boot)
    return scores, scoring.best_pooled_dict(scores)


if __name__ == "__main__":
    paths = importlib.import_module("paths")
    load_validate = importlib.import_module("01_load_validate")

    df_clean = pd.read_csv(paths.table_path("df_clean.csv"))
    valid_groups = [g for g in load_validate.GROUP_ORDER if g in set(df_clean["Group2"])]
    df_eval, excl, valid_eval = restrict_to_evaluable(df_clean, valid_groups, load_validate.GROUP_ORDER)
    df_eval.to_csv(paths.table_path("df_eval.csv"), index=False)
    excl.to_csv(paths.table_path("rows_excluded_not_evaluable.csv"), index=False)
    print(excl.to_string(index=False))

    results_long = run_formula_methods(df_eval, valid_eval)
    scores_per_group, _ = score_formula_methods_per_group(results_long)
    scores_pooled, best_pooled = score_formula_methods_pooled(results_long)
    results_long.to_csv(paths.table_path("formula_results_long.csv"), index=False)
    scores_per_group.to_csv(paths.table_path("formula_scores_per_group.csv"), index=False)
    scores_pooled.to_csv(paths.table_path("formula_scores_pooled.csv"), index=False)

    print("\nPOOLED formula scores (headline summary)")
    print(scores_pooled.sort_values(["Base", "Pearson_r"], ascending=[True, False])
          [["Method", "Base", "Pearson_r", "Pearson_r_lo", "Pearson_r_hi", "CCC", "MAE", "N"]].to_string(index=False))
