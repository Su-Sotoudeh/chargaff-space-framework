# Chargaff composition prediction — Figure 7 pipeline (v2.1, organism-level evaluation)

Formula methods (M1–M5) versus four ML models (RandomForest, GradientBoost, SVR,
NeuralNet) for predicting whole-chromosome composition (A_W, T_W, G_W, C_W),
scored independently for every Group2.

## What changed in v2 and why

### 1. Train/test split is by ORGANISM, never by chromosome (the editor's point)
Chromosomes of one organism share GC content, Z, annotation composition, etc.
Training on human chr1 and testing on chr2 therefore lets the model "see" the
test organism. In v2:

* the split unit is `Organ_ID`; all chromosomes of an organism are either all
  in training or all in testing (`grouped_split.py`, asserts no leakage);
* `StratifiedGroupKFold`, 5 folds, strata = organism-mean GC quartile, so every
  test fold covers low- and high-GC organisms;
* repeated 3 times with different seeds (mean and SD over repeats);
* every row receives exactly one out-of-fold prediction per repeat;
* big organisms (many chromosomes) can no longer dominate: the bootstrap CIs
  also resample organisms, not rows;
* optional: set `SPLIT_GROUP_COL` in `03_ml_methods.py` to a taxonomy column
  (genus/family) for a stricter leave-clade-out test.

Per-group independence is kept: a model trained on Chordata never sees Plasmids.
Groups with < 10 organisms or < 30 rows are skipped for ML only (an
organism-level 5-fold split is meaningless with e.g. 6 organisms); formula
methods are still scored there.

### 2. Headline methods are scored on identical rows (fixed in v2.1)
Rows are kept when the observed composition, the baseline M0 and the
dominant-strand formula M4 are finite; ML is trained/tested on exactly these
rows. **Single-stranded rows (Z = 0 or 100: ssDNA/ssRNA viruses, many plasmids)
are KEPT.** (v2.0 wrongly dropped every row where *any* formula was undefined,
which removed DNA Viruses(+), DNA Viruses(-), RNA Viruses(+), RNA Viruses(-)
and 44% of Plasmids.)
M1, M2, M3, M5 each need one specific strand, so they are undefined where that
strand is empty (M1 at Z = 0, M2 at Z = 100 ...); they are scored on the rows
where they are defined (see the N column) and their skill vs M0 uses those same
rows. The Friedman test drops such a method (printed) rather than dropping whole
groups. ML now uses both annotated strands (`FEATURE_SET = "plus_and_comp"`,
same information as M4); `"plus_only"` gives ML only the +annotation.

### 3. Baseline and better metrics
* **M0_Identity** baseline (composition of the dominant annotated strand,
  normalised, no Chargaff correction = M4 without the correction). `MAE_skill_vs_M0 = 1 - MAE/MAE_M0` and `Pearson_gain_vs_M0`
  tell you whether a method beats doing nothing.
* Pearson r, Lin's CCC, MAE, RMSE, R², plus
  **Within_org_r** (organism mean removed: can the method predict how a
  chromosome differs from its own organism?) and **Between_org_r**.
* 95% CIs by bootstrapping organisms (Pearson r, CCC, MAE, Within_org_r).

### 4. Statistics
* Winner per base = lowest mean rank across groups (not pooled r, which large
  groups dominate).
* Friedman test + Nemenyi critical-difference bar on the CD diagrams. Caveat:
  (group × base) blocks are not fully independent, so treat the p-value as
  descriptive support, not a strict test.

### 5. Leakage audit (`08_leakage_audit.py`)
Same RandomForest scored with random-row KFold versus organism-level folds →
`tables/leakage_audit.csv`. Report this in the Supplement as a direct answer to
the editor. Interpret the drop in Pearson r / rise in MAE; Within_org_r is a
separate diagnostic and is not a leakage indicator.

## How to run

```bash
cd code
python3 run_pipeline.py --n-jobs 40      # reads Dataframes/Cleaned_datafram_whole_organims_2026_07.csv
# or: python3 run_pipeline.py --csv /path/to/Cleaned_datafram_whole_organims_2026_07.csv --n-jobs 40
```

Options: `--n-splits 5 --n-repeats 3 --min-organisms 10 --n-boot 1000
--n-boot-pooled 200 --no-leakage-audit --no-save-long`.
The CSV must contain `Organ_ID`, `Chr_ID`, `Group2`, `Z`, `A_W..C_W`, and the
`annot_*` columns used before. `00_make_synthetic_test_data.py` creates a small
fake file to smoke-test the code (its numbers are NOT results).

Main figure module: set `MAIN_FIGURE_MODULE` in `run_pipeline.py`
(`"Main_figure_2"` default, or `"06_main_figure"`).

## Pipeline steps

1. load + validate (01) 2. common evaluable rows (02) 3. formula + baseline (02)
4. ML out-of-fold CV (03) 5. scoring with CIs (scoring.py)
6. Friedman / CD (07) 7. supplementary heatmaps (05; Pearson r and within-organism r)
8. main figure (Main_figure_2 / 06) 9. manuscript summary table 10. leakage audit (08)

`tables/run_settings.json` records the settings used for each run.
All outputs go to `Figure_7/outputs/tables/` and `Figure_7/outputs/figures/` (git-ignored).

## Runtime — NOT measured on your real data
On the synthetic test (about 1.8k rows, 2 cores) ML took ~163 s. Scoring at
4.6M long-table rows took ~53 s per-group and ~13 s pooled. The real run
(~110k rows, 5×3 folds × 4 models × 4 bases per group, 48 cores) is expected to
take on the order of 5–15+ minutes, but that is an estimate only. The ML stage
prints progress; reduce `--n-repeats` to 1 for a first pass.

## Notes
* SVR subsamples to 5,000 training rows per fit; other models use all rows.
* No hyperparameter tuning; if you add it, nest it inside the training
  organisms only.
* Gray "n/a" cells in heatmaps: correlation is undefined when predictions or
  observations are constant.

## Files
```
paths.py                     output-path resolver
grouped_split.py   (new)     organism-level StratifiedGroupKFold
scoring.py         (new)     fast scoring + organism bootstrap CIs
00_make_synthetic_test_data.py  smoke-test data only
01_load_validate.py          loading + validation (needs Organ_ID, Chr_ID)
02_formula_methods.py        M0 baseline + M1–M5, common rows
03_ml_methods.py             4 ML models, organism-level repeated CV, parallel
04_plot_style.py             style
05_supplementary_figure.py   per-group heatmaps
06_main_figure.py / Main_figure_2.py   main figure
07_stats_cd_diagram.py       Friedman, Nemenyi, CD diagrams
08_leakage_audit.py  (new)   random vs organism-level split
run_pipeline.py              runs everything
```
