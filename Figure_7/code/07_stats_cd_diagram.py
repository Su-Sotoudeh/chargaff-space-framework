"""
07_stats_cd_diagram.py
------------------------
Friedman test + Critical-Difference (Nemenyi) ranking diagrams, for the
formula-methods per-group scores, the ML per-group scores, or both together.

WHAT THIS TESTS:
Given a per-group, per-base Pearson r score for each method, this builds
"blocks" = (Group2, Base) pairs and "treatments" = methods. The Friedman
test asks: across all these blocks, are the rank differences between
methods statistically significant? (Demsar 2006). The Nemenyi critical
difference (CD) is drawn as a bar: two methods whose average ranks differ by
less than CD cannot be called different at alpha = 0.05.

CAVEATS to state in Methods:
  * the 4 bases of one group are not independent blocks (A+T+G+C = 100),
    so the effective number of blocks is smaller than "n blocks"; treat
    the p-value and CD as indicative, not exact;
  * ML scores are organism-level out-of-fold scores (03_ml_methods.py), formula
    scores are computed on the very same rows (02_formula_methods.py).

NEW in this version
  * nemenyi_cd()                    -> literal CD value, drawn in the diagrams
  * pick_winners_by_mean_rank()     -> picks the "best" method per base by mean
                                      rank across groups (not by pooled r, which
                                      is dominated by the biggest groups)
  * make_cd_diagram()               -> same non-colliding label layout as the
                                      combined diagram (old collision logic removed)

The Friedman test requires COMPLETE blocks: every method must have a
score for every (Group2, Base) block used. Blocks with any missing
method (e.g. ML skipped a group with too few organisms) are dropped from the
test automatically, with a count printed.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import friedmanchisquare, rankdata
import importlib

style = importlib.import_module("04_plot_style")
style.set_publication_style()

# Nemenyi q_alpha (alpha = 0.05), k = 2..10 methods (Demsar 2006, Table 5)
_Q_ALPHA_005 = {2: 1.960, 3: 2.343, 4: 2.569, 5: 2.728, 6: 2.850,
                7: 2.949, 8: 3.031, 9: 3.102, 10: 3.164}


def nemenyi_cd(n_methods, n_blocks):
    """Critical difference in average-rank units; None if k is outside 2..10."""
    q = _Q_ALPHA_005.get(int(n_methods))
    if q is None or n_blocks < 1:
        return None
    return q * np.sqrt(n_methods * (n_methods + 1) / (6.0 * n_blocks))


def pick_winners_by_mean_rank(scores_per_group, methods, value_col="Pearson_r"):
    """
    For each base, the winner is the method with the lowest MEAN RANK across
    groups (rank 1 = highest value_col within a group). Only groups where every
    candidate has a score are used (falls back to all groups if none qualify).
    Returns {base: method}.
    """
    d = scores_per_group[scores_per_group["Method"].isin(methods)].dropna(subset=[value_col])
    winners = {}
    for base, sub in d.groupby("Base", observed=True):
        piv = sub.pivot(index="Group2", columns="Method", values=value_col)
        complete = piv.dropna()
        use = complete if len(complete) else piv
        ranks = use.rank(axis=1, ascending=False, method="average")
        winners[base] = str(ranks.mean(axis=0).idxmin())
    return winners


def combine_formula_and_ml_scores(formula_scores_per_group, ml_scores_per_group,
                                    method_col="Method", value_col="Pearson_r",
                                    block_cols=("Group2", "Base")):
    """
    Stacks the formula and ML per-group score tables into one long table
    with all 9 methods (5 formula + 4 ML) so they can be ranked against each
    other within the SAME (Group2, Base) blocks. Blocks where ML skipped the
    group (too few organisms) are incomplete and are excluded automatically.
    """
    keep_cols = list(block_cols) + [method_col, value_col]
    f = formula_scores_per_group[keep_cols].copy()
    m = ml_scores_per_group[keep_cols].copy()
    for df_ in (f, m):
        for c in block_cols + (method_col,):
            df_[c] = df_[c].astype(str)
    return pd.concat([f, m], ignore_index=True)


# methods that are defined for every row (never dropped to rescue blocks)
NEVER_DROP = {"M0_Identity", "M4_Dominant", "RandomForest", "GradientBoost", "SVR", "NeuralNet"}


def run_friedman_and_ranks(scores_per_group, method_col="Method", value_col="Pearson_r",
                            block_cols=("Group2", "Base"), min_block_frac=0.9):
    """
    Returns: friedman_stat, friedman_p, avg_ranks (Series, sorted best-first),
             n_blocks_used, n_blocks_total

    The Friedman test needs COMPLETE blocks. Some formula methods are undefined
    for single-stranded groups (M1 at Z=0, M2 at Z=100 ...). If keeping only
    complete blocks would discard more than (1 - min_block_frac) of all blocks
    (i.e. whole groups such as ssRNA viruses), the method with the most missing
    blocks is dropped from the test instead (printed), until enough blocks
    remain. Set min_block_frac=1.0 to always keep every method.
    """
    scores_per_group = scores_per_group.copy()
    for c in tuple(block_cols) + (method_col,):
        scores_per_group[c] = scores_per_group[c].astype(str)
    df = scores_per_group.dropna(subset=[value_col]).copy()
    methods_with_any_data = sorted(df[method_col].unique())
    methods_all = sorted(scores_per_group[method_col].unique())

    dropped_methods = set(methods_all) - set(methods_with_any_data)
    if dropped_methods:
        print(f"[WARN] Method(s) {sorted(dropped_methods)} had NO valid (non-NaN) "
              f"{value_col} in ANY block and are fully excluded from the Friedman "
              f"test and CD diagram below.")

    methods = methods_with_any_data
    piv = df.pivot_table(index=list(block_cols), columns=method_col, values=value_col, aggfunc="first")
    n_blocks_total = len(piv)
    piv = piv.reindex(columns=methods)
    while piv.shape[1] > 3 and piv.dropna().shape[0] < min_block_frac * n_blocks_total:
        cand = piv.isna().sum()[[c for c in piv.columns if c not in NEVER_DROP]]
        if len(cand) == 0 or cand.max() == 0:
            break          # remaining incomplete blocks are ML-skipped groups: those blocks are excluded
        worst = cand.idxmax()
        print(f"[INFO] Friedman: '{worst}' is undefined in {int(piv[worst].isna().sum())} of "
              f"{n_blocks_total} blocks (single-stranded groups?) -> dropped from this test so those "
              f"groups stay in it.")
        piv = piv.drop(columns=[worst])
    methods = list(piv.columns)
    piv = piv.dropna()
    block_matrix = piv.to_numpy()
    n_blocks_used = len(block_matrix)

    if n_blocks_used < 2 or block_matrix.shape[1] < 3:
        print(f"[WARN] Not enough complete blocks ({n_blocks_used}) or methods "
              f"({block_matrix.shape[1] if len(block_matrix) else 0}) to run Friedman test.")
        return None, None, None, n_blocks_used, n_blocks_total

    print(f"[INFO] Friedman test using methods: {methods} "
          f"({n_blocks_used}/{n_blocks_total} blocks had a complete score for all of them)")

    friedman_stat, friedman_p = friedmanchisquare(
        *[block_matrix[:, i] for i in range(block_matrix.shape[1])]
    )
    rank_matrix = np.array([rankdata(-row, method="average") for row in block_matrix])
    avg_ranks = pd.Series(rank_matrix.mean(axis=0), index=methods).sort_values()
    return friedman_stat, friedman_p, avg_ranks, n_blocks_used, n_blocks_total


def _draw_cd_bar(ax, cd, xmin, y_extent, side, fs, compact):
    """CD bar placed on the less crowded side (side=+1 top, -1 bottom)."""
    if cd is None:
        return
    y_bar = side * (y_extent - (0.45 if compact else 0.3))
    ax.hlines(y_bar, xmin, xmin + cd, color="black", linewidth=3.5, zorder=3)
    ax.vlines([xmin, xmin + cd], y_bar - 0.08, y_bar + 0.08, color="black", linewidth=1.2, zorder=3)
    ax.text(xmin + cd / 2, y_bar - side * 0.12, f"CD = {cd:.2f}", ha="center",
            va="top" if side > 0 else "bottom", fontsize=fs, fontweight="bold")


def make_combined_cd_diagram(avg_ranks, friedman_p, n_blocks_used,
                              formula_methods, ml_methods, title, out_path=None,
                              ax=None, compact=False, show_cd=True):
    """
    Formula methods (blue circles) always ABOVE the rank line, ML methods
    (vermillion triangles) always BELOW. If ax is provided, draws into that
    axes (embedding in a multi-panel figure); otherwise creates a standalone
    figure and saves it when out_path is given. show_cd draws the Nemenyi CD bar.
    """
    n_methods = len(avg_ranks)
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(9, 1.4 + 0.5 * n_methods))
    else:
        fig = ax.figure

    fs_label = 7.0 if compact else 8.5
    fs_tick = 7.5 if compact else 9
    fs_title = 8.5 if compact else 10.5
    marker_size = 7 if compact else 10
    label_template = "{m}\n({r:.2f})" if compact else "{m}\n(rank={r:.2f})"

    xmin, xmax = 1, n_methods
    ax.set_xlim(xmin - 0.4, xmax + 0.4)
    depth_step = 0.78 if compact else 0.6
    formula_color = style.OKABE_ITO["blue"]
    ml_color = style.OKABE_ITO["vermillion"]
    is_formula_arr = np.array([m in formula_methods for m in avg_ranks.index])
    n_formula_present = int(is_formula_arr.sum())
    n_ml_present = n_methods - n_formula_present
    max_depth_per_side = max(n_formula_present, n_ml_present, 1)
    y_extent = max_depth_per_side * depth_step + (1.15 if compact else 0.85)
    ax.set_ylim(-y_extent, y_extent)

    ax.hlines(0, xmin, xmax, color="black", linewidth=1.1, zorder=2)
    ax.set_xticks([])
    tick_drop = 0.22 if compact else 0.16
    for i in range(1, n_methods + 1):
        ax.vlines(i, -tick_drop, 0, color="black", linewidth=1.0, zorder=2)
        ax.text(i, -tick_drop - (0.14 if compact else 0.1), str(i),
                 ha="center", va="top", fontsize=fs_tick)

    rank_vals = avg_ranks.values
    sort_order = np.argsort(rank_vals)
    side = np.array([1 if is_formula_arr[i] else -1 for i in range(n_methods)])
    depth = np.zeros(n_methods, dtype=int)
    side_counts = {1: 0, -1: 0}
    for idx in sort_order:
        s = side[idx]
        side_counts[s] += 1
        depth[idx] = side_counts[s]

    for k, (method, rank_val) in enumerate(avg_ranks.items()):
        is_formula = method in formula_methods
        color = formula_color if is_formula else ml_color
        marker = "o" if is_formula else "^"
        y_text = side[k] * (depth_step * depth[k] + (0.18 if compact else 0.12))
        ax.vlines(rank_val, 0, y_text * 0.85, color="dimgray", linestyle=(0, (2, 1.5)),
                   alpha=0.85, linewidth=1.0, zorder=1)
        ax.plot(rank_val, 0, marker, markersize=marker_size, color=color,
                 markeredgecolor="black", markeredgewidth=0.6, zorder=3)
        va = "bottom" if side[k] > 0 else "top"
        ax.text(rank_val, y_text, label_template.format(m=method, r=rank_val),
                 ha="center", va=va, fontsize=fs_label, fontweight="bold")

    if show_cd:
        cd = nemenyi_cd(n_methods, n_blocks_used)
        _draw_cd_bar(ax, cd, xmin, y_extent, -1 if n_ml_present <= n_formula_present else 1, fs_label, compact)

    formula_handle = plt.Line2D([], [], marker="o", linestyle="", color=formula_color,
                                  markeredgecolor="black", markeredgewidth=0.6,
                                  markersize=marker_size, label="Formula method")
    ml_handle = plt.Line2D([], [], marker="^", linestyle="", color=ml_color,
                             markeredgecolor="black", markeredgewidth=0.6,
                             markersize=marker_size, label="ML method")
    ax.legend(handles=[formula_handle, ml_handle], loc="upper center",
              bbox_to_anchor=(0.5, -0.28 if compact else -0.30),
              ncol=2, fontsize=fs_label, frameon=False)

    ax.set_yticks([])
    ax.set_xlabel("Average rank across all (group, base) blocks  (1 = best)",
                  fontsize=fs_tick, labelpad=10 if compact else 8)
    p_str = f"{friedman_p:.2e}" if friedman_p is not None else "N/A"
    ax.set_title(f"{title}\nFriedman p = {p_str}  |  n blocks = {n_blocks_used}",
                 fontsize=fs_title, pad=14 if compact else 8)
    for i in range(1, n_methods + 1):
        ax.axvline(i, color="gray", alpha=0.15, linewidth=0.7, zorder=0)

    if standalone and out_path is not None:
        style.save_fig(fig, out_path)
    return fig, ax


def make_cd_diagram(avg_ranks, friedman_p, n_blocks_used, title, out_path, cd_value="auto"):
    """
    Single-family ranking diagram. Labels alternate above/below the line and
    every label gets its own row on its side, so none can collide.
    cd_value: "auto" -> Nemenyi CD from (n_methods, n_blocks_used); a number
    -> use it; None -> no CD bar.
    """
    n_methods = len(avg_ranks)
    fig, ax = plt.subplots(figsize=(9, 1.4 + 0.5 * n_methods))
    xmin, xmax = 1, n_methods
    ax.set_xlim(xmin - 0.4, xmax + 0.4)
    depth_step = 0.6
    per_side = int(np.ceil(n_methods / 2))
    y_extent = per_side * depth_step + 1.1
    ax.set_ylim(-y_extent, y_extent)

    ax.hlines(0, xmin, xmax, color="black", linewidth=1.2)
    ax.set_xticks([])
    for i in range(1, n_methods + 1):
        ax.vlines(i, -0.16, 0, color="black", linewidth=1.0)
        ax.text(i, -0.26, str(i), ha="center", va="top", fontsize=9)
        ax.axvline(i, color="gray", alpha=0.15, linewidth=0.7, zorder=0)

    colors = [style.OKABE_ITO[k] for k in
              ["blue", "vermillion", "bluish_green", "orange", "reddish_purple",
               "sky_blue", "yellow", "gray", "cyan"]]
    side_counts = {1: 0, -1: 0}
    for k, (method, rank_val) in enumerate(avg_ranks.items()):
        s = 1 if k % 2 == 0 else -1
        side_counts[s] += 1
        y_text = s * (depth_step * side_counts[s] + 0.12)
        ax.vlines(rank_val, 0, y_text * 0.85, color="dimgray", linestyle=(0, (2, 1.5)), alpha=0.85, linewidth=1.0)
        ax.plot(rank_val, 0, "o", markersize=10, color=colors[k % len(colors)],
                 markeredgecolor="black", zorder=3)
        ax.text(rank_val, y_text, f"{method}\n(rank={rank_val:.2f})", ha="center",
                 va="bottom" if s > 0 else "top", fontsize=8.5, fontweight="bold")

    cd = nemenyi_cd(n_methods, n_blocks_used) if cd_value == "auto" else cd_value
    _draw_cd_bar(ax, cd, xmin, y_extent, 1, 9, False)

    ax.set_yticks([])
    ax.set_xlabel("Average rank across all (group, base) blocks  (1 = best)")
    p_str = f"{friedman_p:.2e}" if friedman_p is not None else "N/A"
    ax.set_title(f"{title}\nFriedman p = {p_str}  |  n blocks = {n_blocks_used}", fontsize=10.5)
    style.save_fig(fig, out_path)
    return fig


if __name__ == "__main__":
    paths = importlib.import_module("paths")
    FORMULA_METHODS = ["M1_PlusOnly", "M2_CompOnly", "M3_Average", "M4_Dominant", "M5_Weighted"]
    ML_METHODS = ["RandomForest", "GradientBoost", "SVR", "NeuralNet"]

    formula_scores_per_group = pd.read_csv(paths.table_path("formula_scores_per_group.csv"))
    ml_scores_per_group = pd.read_csv(paths.table_path("ml_scores_per_group.csv"))

    f_stat, f_p, f_ranks, f_used, f_total = run_friedman_and_ranks(formula_scores_per_group)
    if f_ranks is not None:
        print(f"Formula methods Friedman: stat={f_stat:.3f}, p={f_p:.3e}, blocks used={f_used}/{f_total}")
        print(f_ranks)
        make_cd_diagram(f_ranks, f_p, f_used, "Formula methods: ranking across groups and bases",
                         paths.figure_path("CD_diagram_formula"))

    m_stat, m_p, m_ranks, m_used, m_total = run_friedman_and_ranks(ml_scores_per_group)
    if m_ranks is not None:
        print(f"\nML methods: Friedman stat={m_stat:.3f}, p={m_p:.3e}, blocks used={m_used}/{m_total}")
        print(m_ranks)
        make_cd_diagram(m_ranks, m_p, m_used, "ML methods: ranking across groups and bases",
                         paths.figure_path("CD_diagram_ml"))

    combined_scores = combine_formula_and_ml_scores(formula_scores_per_group, ml_scores_per_group)
    c_stat, c_p, c_ranks, c_used, c_total = run_friedman_and_ranks(combined_scores)
    if c_ranks is not None:
        print(f"\nCOMBINED (formula+ML) Friedman: stat={c_stat:.3f}, p={c_p:.3e}, "
              f"blocks used={c_used}/{c_total}")
        print(c_ranks)
        make_combined_cd_diagram(c_ranks, c_p, c_used, FORMULA_METHODS, ML_METHODS,
                                  "Combined ranking: formula vs ML methods",
                                  out_path=paths.figure_path("CD_diagram_combined"))
    else:
        print("\n[INFO] Combined formula+ML CD diagram skipped (not enough complete blocks).")
