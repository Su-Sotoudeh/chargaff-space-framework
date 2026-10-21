"""
run_pipeline.py
-----------------
Master script (v2: organism-level evaluation).

  1. Load + validate (Organ_ID / Chr_ID required)
  2. Restrict to rows where observed, M0 and M4 are finite -> df_eval (single-stranded rows are kept).
     Formula methods and ML are scored on exactly these rows.
  3. Formula methods M1-M5 + identity baseline M0; scores with organism-bootstrap CIs
  4. ML methods: organism-level, stratified, repeated 5-fold CV, per Group2
  5. ML scores (same scoring code) + skill vs the M0 baseline
  6. Winners per base = lowest mean rank across groups; Friedman + Nemenyi CD diagrams
  7. Supplementary heatmaps (Pearson r, and within-organism r)
  8. Main figure
  9. Manuscript summary table
 10. Leakage audit: random-row split vs organism-level split (optional)

Usage (from the code/ folder):
    python3 run_pipeline.py                         # uses Dataframes/Cleaned_datafram_whole_organims_2026_07.csv
    python3 run_pipeline.py --csv /path/file.csv
    python3 run_pipeline.py --n-repeats 1           # quick run
    python3 run_pipeline.py --n-jobs 40 --no-leakage-audit
"""

import argparse
import importlib
import json
import time

import numpy as np
import pandas as pd

paths = importlib.import_module("paths")
load_validate = importlib.import_module("01_load_validate")
formula_mod = importlib.import_module("02_formula_methods")
ml_mod = importlib.import_module("03_ml_methods")
sup_fig_mod = importlib.import_module("05_supplementary_figure")
stats_mod = importlib.import_module("07_stats_cd_diagram")
scoring = importlib.import_module("scoring")

# Which main-figure layout to draw: "Main_figure_2" (8 panels A-H, newest) or "06_main_figure" (7 panels)
MAIN_FIGURE_MODULE = "Main_figure_2"

BASES = ["A", "T", "G", "C"]
FORMULA_METHODS = formula_mod.FORMULA_METHODS
ML_METHODS = ml_mod.ML_METHODS


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=None, help="path to the genome CSV (default: the master table in Dataframes/, see config.py)")
    ap.add_argument("--n-splits", type=int, default=ml_mod.N_SPLITS)
    ap.add_argument("--n-repeats", type=int, default=ml_mod.N_REPEATS)
    ap.add_argument("--n-jobs", type=int, default=ml_mod.N_JOBS)
    ap.add_argument("--min-organisms", type=int, default=ml_mod.MIN_ORGANISMS)
    ap.add_argument("--n-boot", type=int, default=1000, help="organism-bootstrap resamples for per-group CIs")
    ap.add_argument("--n-boot-pooled", type=int, default=200)
    ap.add_argument("--no-leakage-audit", action="store_true")
    ap.add_argument("--no-save-long", action="store_true", help="do not write the (large) long prediction tables")
    return ap.parse_args()


def audit_oof_coverage(ml_long, df_eval, n_repeats):
    """Every evaluated row must have exactly n_repeats out-of-fold predictions per (method, base)."""
    expected = df_eval["Group2"].value_counts()
    cnt = ml_long.groupby(["Group2", "Method", "Base"], observed=True).size()
    bad = [(k, v, expected[k[0]] * n_repeats) for k, v in cnt.items() if v != expected[k[0]] * n_repeats]
    if bad:
        raise RuntimeError(f"Out-of-fold coverage problem (group, method, base): {bad[:5]}")
    print(f"[OK] Coverage audit: every evaluated row has exactly {n_repeats} out-of-fold prediction(s) "
          f"per method/base, made by models that never saw its organism.")


def matched_baseline_scores(baseline_long, method_long, pooled=False):
    """M0 scored, for every formula method, on exactly the rows where that method
    is defined (M1/M2/M3/M5 are undefined where their strand is empty)."""
    frames = []
    for m in method_long["Method"].astype(str).unique():
        rid = method_long.loc[method_long["Method"].astype(str) == m, "RowID"].unique()
        b = baseline_long[baseline_long["RowID"].isin(rid)].copy()
        b["Method"] = m
        frames.append(b)
    bl = pd.concat(frames, ignore_index=True)
    return scoring.score_pooled(bl, n_boot=0) if pooled else scoring.score_per_group(bl, n_boot=0)


def main():
    args = parse_args()
    main_fig_mod = importlib.import_module(MAIN_FIGURE_MODULE)
    t_start = time.time()

    print("\n" + "#" * 70 + "\nSTEP 1: LOAD + VALIDATE\n" + "#" * 70)
    csv_path = args.csv or load_validate.CSV_PATH
    df_clean, valid_groups = load_validate.load_and_validate(csv_path=csv_path)
    df_clean.to_csv(paths.table_path("df_clean.csv"), index=False)
    print(f"[TIME] Step 1 done at {time.time()-t_start:.1f}s")

    print("\n" + "#" * 70 + "\nSTEP 2: COMMON EVALUATION ROWS (M0, M4 and ML scored on the same rows)\n" + "#" * 70)
    df_eval, excl, valid_eval = formula_mod.restrict_to_evaluable(df_clean, valid_groups, load_validate.GROUP_ORDER)
    df_eval.to_csv(paths.table_path("df_eval.csv"), index=False)
    excl.to_csv(paths.table_path("rows_excluded_not_evaluable.csv"), index=False)
    print(excl.to_string(index=False))
    print(f"[INFO] {len(df_eval)} rows / {df_eval['UnitCode'].nunique()} organism-units kept (headline methods M0, M4 and all ML are scored on exactly these rows).")

    print("\n" + "#" * 70 + "\nSTEP 3: FORMULA METHODS + BASELINE (independent per group)\n" + "#" * 70)
    formula_long = formula_mod.run_formula_methods(df_eval, valid_eval)
    baseline_long = formula_mod.run_baseline_identity(df_eval, valid_eval)
    formula_scores_per_group, _ = formula_mod.score_formula_methods_per_group(formula_long, n_boot=args.n_boot)
    formula_scores_pooled, _ = formula_mod.score_formula_methods_pooled(formula_long, n_boot=args.n_boot_pooled)
    base_per_group = scoring.score_per_group(baseline_long, n_boot=0)
    base_pooled = scoring.score_pooled(baseline_long, n_boot=0)
    formula_scores_per_group = scoring.add_baseline_skill(
        formula_scores_per_group, matched_baseline_scores(baseline_long, formula_long), ["Group2", "Method", "Base"])
    formula_scores_pooled = scoring.add_baseline_skill(
        formula_scores_pooled, matched_baseline_scores(baseline_long, formula_long, pooled=True), ["Method", "Base"])
    best_formula_per_base = stats_mod.pick_winners_by_mean_rank(formula_scores_per_group, FORMULA_METHODS)
    print("\nPOOLED formula scores (headline only):")
    print(formula_scores_pooled.sort_values(["Base", "Pearson_r"], ascending=[True, False])
          [["Method", "Base", "Pearson_r", "Pearson_r_lo", "Pearson_r_hi", "CCC", "MAE", "MAE_skill_vs_M0", "N"]]
          .round(3).to_string(index=False))
    print("Best formula per base (lowest mean rank across groups):", best_formula_per_base)
    print(f"[TIME] Step 3 done at {time.time()-t_start:.1f}s")

    print("\n" + "#" * 70 + "\nSTEP 4: ML METHODS (organism-level repeated CV, independent per group)\n" + "#" * 70)
    ml_long, skipped_groups = ml_mod.run_ml_methods_per_group(
        df_eval, valid_eval, n_splits=args.n_splits, n_repeats=args.n_repeats,
        min_organisms=args.min_organisms, n_jobs=args.n_jobs)
    ran_groups = sorted(set(ml_long["Group2"].astype(str)))
    audit_oof_coverage(ml_long, df_eval[df_eval["Group2"].isin(ran_groups)], args.n_repeats)

    print("\n" + "#" * 70 + "\nSTEP 5: ML SCORES\n" + "#" * 70)
    ml_scores_per_group, _ = ml_mod.score_ml_methods_per_group(ml_long, n_boot=args.n_boot)
    ml_scores_pooled, _ = ml_mod.score_ml_methods_pooled(ml_long, n_boot=args.n_boot_pooled)
    ml_scores_per_group = scoring.add_baseline_skill(ml_scores_per_group, base_per_group, ["Group2", "Base"])
    ml_scores_pooled = scoring.add_baseline_skill(ml_scores_pooled, base_pooled, ["Base"])
    best_ml_per_base = stats_mod.pick_winners_by_mean_rank(ml_scores_per_group, ML_METHODS)
    print(ml_scores_pooled.sort_values(["Base", "Pearson_r"], ascending=[True, False])
          [["Method", "Base", "Pearson_r", "Pearson_r_lo", "Pearson_r_hi", "CCC", "MAE", "MAE_skill_vs_M0", "N"]]
          .round(3).to_string(index=False))
    print("Best ML per base (lowest mean rank across groups):", best_ml_per_base)

    # ---- save tables ----
    formula_scores_per_group.to_csv(paths.table_path("formula_scores_per_group.csv"), index=False)
    formula_scores_pooled.to_csv(paths.table_path("formula_scores_pooled.csv"), index=False)
    ml_scores_per_group.to_csv(paths.table_path("ml_scores_per_group.csv"), index=False)
    ml_scores_pooled.to_csv(paths.table_path("ml_scores_pooled.csv"), index=False)
    base_per_group.to_csv(paths.table_path("baseline_M0_scores_per_group.csv"), index=False)
    if not args.no_save_long:
        formula_long.to_csv(paths.table_path("formula_results_long.csv"), index=False)
        baseline_long.to_csv(paths.table_path("baseline_M0_results_long.csv"), index=False)
        ml_long.to_csv(paths.table_path("ml_results_long.csv"), index=False)
    with open(paths.table_path("run_settings.json"), "w") as f:
        json.dump({"csv": csv_path, "split_unit": ml_mod.SPLIT_GROUP_COL, "n_splits": args.n_splits,
                   "n_repeats": args.n_repeats, "min_organisms": args.min_organisms,
                   "min_rows": ml_mod.MIN_SAMPLES, "feature_set": ml_mod.FEATURE_SET,
                   "random_state": ml_mod.RANDOM_STATE, "n_boot": args.n_boot,
                   "svr_max_train_rows": ml_mod.SVR_MAX_TRAIN_ROWS_PER_GROUP,
                   "groups_skipped_for_ml": [list(map(str, s)) for s in skipped_groups],
                   "rows_evaluated": int(len(df_eval))}, f, indent=2)
    print(f"[TIME] Step 5 done at {time.time()-t_start:.1f}s")

    print("\n" + "#" * 70 + "\nSTEP 6: FRIEDMAN TEST + NEMENYI CD DIAGRAMS\n" + "#" * 70)
    f_stat, f_p, f_ranks, f_used, _ = stats_mod.run_friedman_and_ranks(formula_scores_per_group)
    if f_ranks is not None:
        print(f"Formula methods: Friedman stat={f_stat:.3f}, p={f_p:.3e}, blocks used={f_used}")
        print(f_ranks)
        stats_mod.make_cd_diagram(f_ranks, f_p, f_used, "Formula methods: ranking across groups and bases",
                                  paths.figure_path("CD_diagram_formula"))
    m_stat, m_p, m_ranks, m_used, _ = stats_mod.run_friedman_and_ranks(ml_scores_per_group)
    if m_ranks is not None:
        print(f"\nML methods: Friedman stat={m_stat:.3f}, p={m_p:.3e}, blocks used={m_used}")
        print(m_ranks)
        stats_mod.make_cd_diagram(m_ranks, m_p, m_used, "ML methods: ranking across groups and bases",
                                  paths.figure_path("CD_diagram_ml"))
    combined = stats_mod.combine_formula_and_ml_scores(formula_scores_per_group, ml_scores_per_group)
    c_stat, c_p, c_ranks, c_used, _ = stats_mod.run_friedman_and_ranks(combined)
    if c_ranks is not None:
        print(f"\nCOMBINED: Friedman stat={c_stat:.3f}, p={c_p:.3e}, blocks used={c_used}, "
              f"Nemenyi CD={stats_mod.nemenyi_cd(len(c_ranks), c_used):.2f}")
        print(c_ranks)
        stats_mod.make_combined_cd_diagram(c_ranks, c_p, c_used, FORMULA_METHODS, ML_METHODS,
                                           "Combined ranking: formula vs ML methods",
                                           out_path=paths.figure_path("CD_diagram_combined"))

    print("\n" + "#" * 70 + "\nSTEP 7: SUPPLEMENTARY FIGURES\n" + "#" * 70)
    for col, lab, tag in (("Pearson_r", "Pearson r", "pearson"), ("Within_org_r", "within-organism r", "within_org_r")):
        sup_fig_mod.make_group_by_method_heatmap_figure(
            formula_scores_per_group, FORMULA_METHODS,
            f"Supplementary Figure S1 ({lab}): formula methods, each group scored independently",
            paths.figure_path(f"SupFig_formula_per_group_{tag}"), value_col=col)
        sup_fig_mod.make_group_by_method_heatmap_figure(
            ml_scores_per_group, ML_METHODS,
            f"Supplementary Figure S2 ({lab}): ML methods on held-out organisms (organism-level CV)",
            paths.figure_path(f"SupFig_ml_per_group_{tag}"), value_col=col)

    print("\n" + "#" * 70 + f"\nSTEP 8: MAIN FIGURE ({MAIN_FIGURE_MODULE})\n" + "#" * 70)
    main_fig_mod.make_main_figure(
        formula_long, ml_long[ml_long["Split"] == "test"],
        formula_scores_per_group, ml_scores_per_group,
        best_formula_per_base, best_ml_per_base,
        paths.figure_path("MainFig_formula_vs_ML"))

    print("\n" + "#" * 70 + "\nSTEP 9: MANUSCRIPT SUMMARY TABLE\n" + "#" * 70)
    rows = []
    for base in BASES:
        fm, mm = best_formula_per_base[base], best_ml_per_base[base]
        fr = formula_scores_pooled[(formula_scores_pooled["Method"] == fm) & (formula_scores_pooled["Base"] == base)].iloc[0]
        mr = ml_scores_pooled[(ml_scores_pooled["Method"] == mm) & (ml_scores_pooled["Base"] == base)].iloc[0]
        rows.append({
            "Base": base,
            "Best_Formula_Method": fm,
            "Formula_Pearson_r_pooled": round(fr["Pearson_r"], 3),
            "Formula_r_95CI": f"{fr['Pearson_r_lo']:.3f}-{fr['Pearson_r_hi']:.3f}",
            "Formula_CCC": round(fr["CCC"], 3), "Formula_MAE": round(fr["MAE"], 3),
            "Formula_within_org_r": round(fr["Within_org_r"], 3),
            "Formula_MAE_skill_vs_M0": round(fr["MAE_skill_vs_M0"], 3),
            "Best_ML_Method": mm,
            "ML_Pearson_r_pooled": round(mr["Pearson_r"], 3),
            "ML_r_95CI": f"{mr['Pearson_r_lo']:.3f}-{mr['Pearson_r_hi']:.3f}",
            "ML_CCC": round(mr["CCC"], 3), "ML_MAE": round(mr["MAE"], 3),
            "ML_within_org_r": round(mr["Within_org_r"], 3),
            "ML_MAE_skill_vs_M0": round(mr["MAE_skill_vs_M0"], 3),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(paths.table_path("manuscript_summary_table.csv"), index=False)
    print(summary.T.to_string(header=False))

    if not args.no_leakage_audit:
        print("\n" + "#" * 70 + "\nSTEP 10: LEAKAGE AUDIT (random-row split vs organism-level split)\n" + "#" * 70)
        audit = importlib.import_module("08_leakage_audit")
        tab = audit.run_leakage_audit(df_eval, valid_eval, n_splits=args.n_splits, n_jobs=args.n_jobs,
                                      min_organisms=args.min_organisms)
        if len(tab):
            tab.to_csv(paths.table_path("leakage_audit.csv"), index=False)
            print(audit.summarize(tab).to_string())

    print(f"\n[TIME] Total pipeline time: {time.time()-t_start:.1f}s")
    print(f"Outputs: {paths.FIGURES_DIR}  |  {paths.TABLES_DIR}")


if __name__ == "__main__":
    main()
