"""
05_supplementary_figure.py
----------------------------
Supplementary Figure: full per-group model-selection evidence.

For each of formula methods and ML methods, produces a 4-panel figure
(one panel per base: A, T, G, C). Each panel is a Group2 x Method
heatmap of Pearson r (or any other column of the score table, e.g. Within_org_r), computed entirely from the PER-GROUP scoring
table -- i.e. each cell reflects only that group's own predictions vs
its own observed values, with no influence from any other group.

This is the "we tried everything honestly, per taxonomic/genome-type
group" evidence panel: it lets you (and reviewers) see directly whether
a method's apparent advantage is consistent across all 15 groups or
driven by just a few.

Two supplementary figures are produced:
    SupFig_formula_per_group_pearson  (5 methods x 15 groups x 4 bases)
    SupFig_ml_per_group_pearson        (4 methods x up-to-15 groups x 4 bases)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import importlib

style = importlib.import_module("04_plot_style")
style.set_publication_style()

BASES = ["A", "T", "G", "C"]

GROUP_ORDER = [
    "Chordata", "Arthropoda", "Plantae", "Fungi", "Archaea", "Bacteria",
    "Chloroplast", "Mitochondria", "Plasmids",
    "DNA Viruses(+/-)", "DNA Viruses(+)", "DNA Viruses(-)",
    "RNA Viruses(+/-)", "RNA Viruses(+)", "RNA Viruses(-)"
]


def make_group_by_method_heatmap_figure(scores_per_group, method_order, title, out_path,
                                         value_col="Pearson_r", vmin=-1, vmax=1, cmap="RdYlGn"):
    """
    scores_per_group: DataFrame with columns [Group2, Method, Base, value_col, ...]
    method_order: list of method names, controls column order
    """
    fig, axes = plt.subplots(2, 2, figsize=(13, 12))
    axes = axes.flatten()

    groups_present = [g for g in GROUP_ORDER if g in scores_per_group["Group2"].unique()]

    for i, base in enumerate(BASES):
        ax = axes[i]
        sub = scores_per_group[scores_per_group["Base"] == base]
        pivot = sub.pivot(index="Group2", columns="Method", values=value_col)
        pivot = pivot.reindex(index=groups_present, columns=method_order)

        im = ax.imshow(pivot.values, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(method_order)))
        ax.set_xticklabels(method_order, rotation=30, ha="right", fontsize=8)
        ax.set_yticks(range(len(groups_present)))
        ax.set_yticklabels(groups_present, fontsize=8)
        ax.set_title(f"Base {base}", fontsize=10.5, fontweight="bold")

        for r in range(pivot.shape[0]):
            for c in range(pivot.shape[1]):
                val = pivot.values[r, c]
                if not np.isnan(val):
                    ax.text(c, r, f"{val:.2f}", ha="center", va="center", fontsize=6.5)
                else:
                    # NaN means the model collapsed to a constant prediction
                    # (zero-variance output -> Pearson r undefined) or the
                    # group/method combination had too few rows to score.
                    # Mark explicitly rather than leaving a silent blank.
                    ax.add_patch(plt.Rectangle((c - 0.5, r - 0.5), 1, 1,
                                                facecolor="lightgray", edgecolor="white"))
                    ax.text(c, r, "n/a", ha="center", va="center", fontsize=6, color="gray", style="italic")

        fig.colorbar(im, ax=ax, shrink=0.75, label=value_col.replace("_", " "))
        style.add_panel_label(ax, ["A", "B", "C", "D"][i], x=-0.35, y=1.05)

    n_nan = scores_per_group[value_col].isna().sum()
    if n_nan > 0:
        print(f"[WARN] {n_nan} (Group, Method, Base) combinations had undefined "
              f"{value_col} -- shown as 'n/a' in the heatmap. Causes: predictions with zero "
              f"variance (too little training data), too few rows, or -- for Within_org_r -- "
              f"groups in which no organism has 2+ chromosomes (e.g. organelles, viruses).")

    fig.suptitle(title, fontsize=12, y=1.01)
    fig.text(0.5, -0.01,
              "'n/a' = undefined (constant predictions, too few rows, or no organism with 2+ "
              "chromosomes for within-organism r)",
              ha="center", fontsize=7.5, style="italic", color="gray")
    fig.tight_layout()
    style.save_fig(fig, out_path)
    return fig


if __name__ == "__main__":
    paths = importlib.import_module("paths")

    formula_scores_per_group = pd.read_csv(paths.table_path("formula_scores_per_group.csv"))
    ml_scores_per_group = pd.read_csv(paths.table_path("ml_scores_per_group.csv"))

    FORMULA_METHODS = ["M1_PlusOnly", "M2_CompOnly", "M3_Average", "M4_Dominant", "M5_Weighted"]
    ML_METHODS = ["RandomForest", "GradientBoost", "SVR", "NeuralNet"]

    make_group_by_method_heatmap_figure(
        formula_scores_per_group, FORMULA_METHODS,
        "Supplementary Figure S1: Formula methods, per-group Pearson r (each group scored independently)",
        paths.figure_path("SupFig_formula_per_group_pearson")
    )

    make_group_by_method_heatmap_figure(
        ml_scores_per_group, ML_METHODS,
        "Supplementary Figure S2: ML methods, per-group Pearson r on held-out ORGANISMS "
        "(organism-level cross-validation; each group trained and scored independently)",
        paths.figure_path("SupFig_ml_per_group_pearson")
    )
