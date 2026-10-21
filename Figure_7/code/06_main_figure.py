"""
06_main_figure.py
-------------------
MAIN FIGURE: Best formula method vs best ML method, across all 4 bases
and all 15 genome groups.

The "best" method per base is chosen by mean rank across groups (see 07_stats_cd_diagram.pick_winners_by_mean_rank)
(one number per method, easy to state in text: "M3_Average was the
best-performing formula method for base A"). But the per-group panel
(E) and the scatter panels (A-D) are built from predictions that were
generated independently per group (formula) or trained independently
per group (ML) -- so no group's own fit quality is contaminated by
another group's data. The pooling only happens at the very end, for
the purpose of choosing one method name and drawing one combined
scatter plot; it never feeds back into how any individual prediction
was produced.

Layout (6 panels):
  A-D: Predicted vs Observed scatter, one per base (A,T,G,C).
       Best formula (blue) and best ML (vermillion) overlaid,
       identity line shown, Pearson r / MAE annotated for each
       (computed on the same pooled-for-display data as the scatter).
  E:   Grouped bar chart -- per-Group2 Pearson r for the winning
       formula vs winning ML method, using the PER-GROUP scores table
       (i.e. the score Su already trusts -- never recomputed by
       pooling here).
  F:   Compact summary heatmap -- winner method's pooled Pearson r per
       base, formula row vs ML row (ties back to the supplementary
       per-group heatmaps).
  G:   Combined critical-difference diagram (built in 07_stats_cd_diagram.py).

===============================================================================
HOW TO ADJUST THIS FIGURE
===============================================================================
Every position, font size, legend location, and spacing value lives in the
CFG dictionary directly below. Find the panel you want to change
(panels_abcd / panel_e / panel_f / text_panel / panel_g / figure), change
the value, re-run. You should never need to edit the plotting code beneath
CFG to move a legend or resize something.

NOTE on Panel G's legend: the critical-difference diagram in Panel G is
drawn by stats_mod.make_combined_cd_diagram(), which lives in
07_stats_cd_diagram.py -- NOT in this file. The panel_g section below only
controls this panel's letter-label position. To move the CD diagram's own
legend, that function needs the same kind of CFG treatment -- share that
file and it can be added there too.
===============================================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error
import importlib

style = importlib.import_module("04_plot_style")
style.set_publication_style()
stats_mod = importlib.import_module("07_stats_cd_diagram")

BASES = ["A", "T", "G", "C"]

GROUP_ORDER = [
    "Chordata", "Arthropoda", "Plantae", "Fungi", "Archaea", "Bacteria",
    "Chloroplast", "Mitochondria", "Plasmids",
    "DNA Viruses(+/-)", "DNA Viruses(+)", "DNA Viruses(-)",
    "RNA Viruses(+/-)", "RNA Viruses(+)", "RNA Viruses(-)"
]

# Biological category color-coding for Panel E's x-axis group labels.
# Organelles = membrane-bound, formerly-endosymbiotic cellular structures
# only (Chloroplast, Mitochondria). Plasmids are mobile genetic elements,
# not organelles, so they're grouped with viruses under "non-cellular."
GROUP_CATEGORY = {
    "Chordata": "Eukaryote", "Arthropoda": "Eukaryote", "Plantae": "Eukaryote", "Fungi": "Eukaryote",
    "Archaea": "Prokaryote", "Bacteria": "Prokaryote",
    "Chloroplast": "Organelle", "Mitochondria": "Organelle",
    "Plasmids": "Non-cellular",
    "DNA Viruses(+/-)": "Non-cellular", "DNA Viruses(+)": "Non-cellular", "DNA Viruses(-)": "Non-cellular",
    "RNA Viruses(+/-)": "Non-cellular", "RNA Viruses(+)": "Non-cellular", "RNA Viruses(-)": "Non-cellular",
}

CATEGORY_COLORS = {
    "Eukaryote": "#1B7837",      # green
    "Prokaryote": "#2166AC",     # blue
    "Organelle": "#B35806",      # brown/orange
    "Non-cellular": "#762A83",   # purple
}


# =============================================================================
# CFG -- every layout / legend / font knob for the figure, grouped by panel.
# Legend positioning uses matplotlib's normal system:
#   loc          = anchor point on the LEGEND BOX itself (e.g. "lower left"
#                  means the legend's lower-left corner sits at bbox_to_anchor)
#   bbox_to_anchor = (x, y) in AXES-FRACTION coordinates for that panel.
#                  (0,0) = bottom-left of the panel, (1,1) = top-right.
#                  y > 1.0 pushes the legend ABOVE the panel; y < 0 pushes it
#                  below. This is the value to raise if a legend sits too low.
# =============================================================================
CFG = {
    "figure": {
        "figsize": (13, 14.5),
        "height_ratios": [1.1, 1.1, 1.65],          # row heights: [A-D row, E row, F/text/G row]
        "hspace": 0.6,
        "wspace": 0.4,
        "top": 0.93,      # top edge of the A-D/E/F/G grid, as a fraction of figure height.
        # Raise toward 1.0 to shrink the blank space above Panel A-D; lower it
        # (e.g. 0.85) to add more room back if the suptitle starts overlapping.
        "bottom_width_ratios": [1.15, 0.95, 2.0],   # F : text : G column widths
        "bottom_wspace": 0.3,
        "suptitle_fontsize": 12,
        "suptitle_y": 1.01,   # if the title still floats too high above Panel A-D after
        # raising "top" above, bring this down too (e.g. 0.97).
    },

    "panels_abcd": {
        "point_size": 9,
        "alpha": 0.55,
        "annotation_fontsize": 6.8,
        "legend_loc": "lower right",
        "legend_bbox": (1.25,0),   # optional (x, y) axes-fraction anchor, e.g. (0, 1.5) puts it
        # well above panel A. None = plain "loc" positioning inside the panel.
        "legend_fontsize": 6.5,
        "title_fontsize": 10,
        "panel_label_x": -0.25,
        "panel_label_y": 1.18,
    },

    "panel_e": {
        "bar_width": 0.35,
        "value_label_fontsize": 6.2,
        "value_label_offset": 0.02,
        "xtick_fontsize": 7.5,
        "ylim": (0.2, 1.15),
        "title_fontsize": 10,
        "title_pad": 38,
        "panel_label_x": -0.05,
        "panel_label_y": 1.28,

        # method-color legend (the blue/vermillion "best formula" / "best ML" key)
        "method_legend_loc": "upper left",
        "method_legend_bbox": (0,1.1),   # optional (x, y) axes-fraction anchor, same idea as above.
        "method_legend_fontsize": 7.5,

        # category-color legend (Eukaryote / Prokaryote / Organelle / Non-cellular).
        # This was previously wrapped in triple-quotes ('''...'''), which turns
        # it into an unused string literal instead of a legend call -- it was
        # never actually being drawn. That bug is fixed below; use
        # show_category_legend to toggle it and the two values under it to
        # move it.
        "show_category_legend": False,
        "category_legend_loc": "lower left",       # anchor point on the legend box
        "category_legend_bbox": (0.0, 1.16),        # (x, y) axes-fraction -- raise y to push it higher
        "category_legend_fontsize": 7.5,
        "category_legend_ncol": 4,
    },

    "panel_f": {
        "width_ratios": [5.5, 0.35],   # heatmap : colorbar
        "wspace": 0.15,
        "cmap": "RdYlGn",
        "vmin": -1,
        "vmax": 1,
        "cell_fontsize": 11,
        "tick_fontsize": 11,
        "cbar_label_fontsize": 8.5,
        "cbar_tick_fontsize": 8,
        "title_fontsize": 10.5,
        "panel_label_x": -0.42,
        "panel_label_y": 1.32,
    },

    "text_panel": {
        "fontsize": 10,
        "linespacing": 1.4,
        "y_start": 0.92,
    },

    "panel_g": {
        "panel_label_x": -0.12,
        "panel_label_y": 1.15,
        # Repositions the CD diagram's legend (built inside
        # stats_mod.make_combined_cd_diagram(), which lives in
        # 07_stats_cd_diagram.py and isn't edited here). We pull the
        # already-drawn legend back out and redraw it at this position,
        # so this works no matter how that function built it.
        # loc = which corner of the LEGEND BOX sits at legend_bbox.
        # legend_bbox = (x, y) in axes-fraction coords, (1,1) = top-right
        # corner of Panel G itself.
        "legend_loc": "upper left",
        "legend_bbox": (0, 1.0),
        "legend_fontsize": 7,
    },
}


def _safe_pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    r, _ = pearsonr(x, y)
    return r


def make_main_figure(formula_long, ml_long_test, formula_scores_per_group,
                     ml_scores_per_group, best_formula_per_base,
                     best_ml_per_base, out_path,
                     formula_methods=None, ml_methods=None, cfg=CFG):
    """
    formula_long              : tidy long df, ALL groups' formula predictions
                                  (each row already produced independently
                                  per group -- this function only pools
                                  them for the display scatter)
    ml_long_test               : tidy long df, ALL groups' ML test
                                  predictions (each group trained
                                  independently -- same note as above)
    formula_scores_per_group   : PER-GROUP formula scores table (the
                                  trusted, non-pooled scoring)
    ml_scores_per_group        : PER-GROUP ML scores table (same)
    best_formula_per_base / best_ml_per_base : dict {base: method_name},
                                  chosen by mean rank across groups (see 07_stats_cd_diagram.pick_winners_by_mean_rank)
    formula_methods / ml_methods : lists of method names for the combined
                                  CD-diagram panel (defaults to the
                                  standard 5 formula + 4 ML method names
                                  if not provided)
    cfg                        : the CFG dict above; pass a modified copy
                                  to override any layout value without
                                  editing this function
    """
    if formula_methods is None:
        formula_methods = ["M1_PlusOnly", "M2_CompOnly", "M3_Average", "M4_Dominant", "M5_Weighted"]
    if ml_methods is None:
        ml_methods = ["RandomForest", "GradientBoost", "SVR", "NeuralNet"]

    # ML long table now holds several CV repeats (organism-level, out-of-fold).
    # For the pooled DISPLAY scatter use ONE repeat, so each chromosome appears once.
    if "Repeat" in ml_long_test.columns:
        ml_long_test = ml_long_test[ml_long_test["Repeat"] == ml_long_test["Repeat"].min()]

    fig_cfg = cfg["figure"]
    abcd_cfg = cfg["panels_abcd"]
    e_cfg = cfg["panel_e"]
    f_cfg = cfg["panel_f"]
    txt_cfg = cfg["text_panel"]
    g_cfg = cfg["panel_g"]

    fig = plt.figure(figsize=fig_cfg["figsize"])
    gs = gridspec.GridSpec(3, 4, height_ratios=fig_cfg["height_ratios"],
                           hspace=fig_cfg["hspace"], wspace=fig_cfg["wspace"], top=fig_cfg["top"])
    # Bottom row gets its own width allocation via a nested GridSpec.
    # Panel F is a 2x4 heatmap -- inherently wide-and-short, so it gets a
    # narrower column than before (it was being stretched tall and thin,
    # the opposite of its natural shape). The winning-method text sits
    # immediately to F's right. Panel G (the CD diagram, needs the most
    # horizontal room for 9 labeled points spread across the full rank
    # range) gets the rest.
    gs_bottom = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=gs[2, :],
        width_ratios=fig_cfg["bottom_width_ratios"], wspace=fig_cfg["bottom_wspace"]
    )

    panel_labels = iter(["A", "B", "C", "D"])

    # ---- Panels A-D: scatter per base (pooled FOR DISPLAY ONLY) ----
    for j, base in enumerate(BASES):
        ax = fig.add_subplot(gs[0, j])

        f_method = best_formula_per_base[base]
        m_method = best_ml_per_base[base]

        f_sub = formula_long[(formula_long["Method"] == f_method) & (formula_long["Base"] == base)]
        m_sub = ml_long_test[(ml_long_test["Method"] == m_method) & (ml_long_test["Base"] == base)]

        ax.scatter(f_sub["Observed"], f_sub["Predicted"], s=abcd_cfg["point_size"], alpha=abcd_cfg["alpha"],
                   color=style.OKABE_ITO["blue"], label=f"Formula: {f_method}",
                   linewidths=0, rasterized=True)
        ax.scatter(m_sub["Observed"], m_sub["Predicted"], s=abcd_cfg["point_size"], alpha=abcd_cfg["alpha"],
                   color=style.OKABE_ITO["vermillion"], label=f"ML: {m_method}",
                   linewidths=0, rasterized=True)

        all_vals = pd.concat([f_sub["Observed"], f_sub["Predicted"],
                              m_sub["Observed"], m_sub["Predicted"]])
        lims = [all_vals.min() - 1, all_vals.max() + 1]
        ax.plot(lims, lims, "--", color="gray", linewidth=1, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)

        r_f = _safe_pearson(f_sub["Observed"], f_sub["Predicted"])
        mae_f = mean_absolute_error(f_sub["Observed"], f_sub["Predicted"])
        r_m = _safe_pearson(m_sub["Observed"], m_sub["Predicted"])
        mae_m = mean_absolute_error(m_sub["Observed"], m_sub["Predicted"])

        ax.text(0.04, 0.97,
                f"Formula r={r_f:.2f}\nMAE={mae_f:.2f}\n\nML r={r_m:.2f}\nMAE={mae_m:.2f}",
                transform=ax.transAxes, fontsize=abcd_cfg["annotation_fontsize"], va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, edgecolor="lightgray"))

        ax.set_xlabel(f"Observed {base} (%)")
        ax.set_ylabel(f"Predicted {base} (%)")
        ax.set_title(f"Base {base}", fontsize=abcd_cfg["title_fontsize"], fontweight="bold")
        if j == 0:
            ax.legend(loc=abcd_cfg["legend_loc"], bbox_to_anchor=abcd_cfg["legend_bbox"],
                      fontsize=abcd_cfg["legend_fontsize"], frameon=False)
        style.add_panel_label(ax, next(panel_labels), x=abcd_cfg["panel_label_x"], y=abcd_cfg["panel_label_y"])

    # ---- Panel E: grouped bars, PER-GROUP Pearson r (trusted table, not pooled) ----
    ax_e = fig.add_subplot(gs[1, :])

    f_rows = []
    m_rows = []
    for base in BASES:
        f_method = best_formula_per_base[base]
        m_method = best_ml_per_base[base]
        f_sub = formula_scores_per_group[
            (formula_scores_per_group["Method"] == f_method) & (formula_scores_per_group["Base"] == base)
            ][["Group2", "Pearson_r"]]
        m_sub = ml_scores_per_group[
            (ml_scores_per_group["Method"] == m_method) & (ml_scores_per_group["Base"] == base)
            ][["Group2", "Pearson_r"]]
        f_rows.append(f_sub)
        m_rows.append(m_sub)

    f_all = pd.concat(f_rows)
    m_all = pd.concat(m_rows)
    # average Pearson r across the 4 bases' winning method, per group
    f_avg = f_all.groupby("Group2")["Pearson_r"].mean()
    m_avg = m_all.groupby("Group2")["Pearson_r"].mean()

    groups_present = [g for g in GROUP_ORDER if g in set(f_avg.index) | set(m_avg.index)]
    x = np.arange(len(groups_present))
    width = e_cfg["bar_width"]

    formula_vals = [f_avg.get(g, np.nan) for g in groups_present]
    ml_vals = [m_avg.get(g, np.nan) for g in groups_present]

    bars_f = ax_e.bar(x - width/2, formula_vals, width, label="Best formula (avg over A,T,G,C)",
                      color=style.OKABE_ITO["blue"])
    bars_m = ax_e.bar(x + width/2, ml_vals, width, label="Best ML (avg over A,T,G,C)",
                      color=style.OKABE_ITO["vermillion"])
    ax_e.axhline(0, color="black", linewidth=0.6)

    # Vertical value labels on top of (or below, if negative) each bar
    for bar, val in zip(bars_f, formula_vals):
        if np.isnan(val):
            continue
        va = "bottom" if val >= 0 else "top"
        offset = e_cfg["value_label_offset"] if val >= 0 else -e_cfg["value_label_offset"]
        ax_e.text(bar.get_x() + bar.get_width() / 2, val + offset, f"{val:.2f}",
                  rotation=90, ha="center", va=va, fontsize=e_cfg["value_label_fontsize"])
    for bar, val in zip(bars_m, ml_vals):
        if np.isnan(val):
            continue
        va = "bottom" if val >= 0 else "top"
        offset = e_cfg["value_label_offset"] if val >= 0 else -e_cfg["value_label_offset"]
        ax_e.text(bar.get_x() + bar.get_width() / 2, val + offset, f"{val:.2f}",
                  rotation=90, ha="center", va=va, fontsize=e_cfg["value_label_fontsize"])

    ax_e.set_xticks(x)
    ax_e.set_xticklabels(groups_present, rotation=45, ha="right", fontsize=e_cfg["xtick_fontsize"])
    # Color-code each group's x-tick label by biological category so
    # eukaryotes / prokaryotes / organelles / non-cellular groups (plasmids
    # + viruses) are visually distinguishable at a glance.
    for tick_label, g in zip(ax_e.get_xticklabels(), groups_present):
        cat = GROUP_CATEGORY.get(g, "Eukaryote")
        tick_label.set_color(CATEGORY_COLORS[cat])
        tick_label.set_fontweight("bold")

    ax_e.set_ylabel("Mean Pearson r (across A,T,G,C)")
    ax_e.set_ylim(e_cfg["ylim"])
    ax_e.set_title("Per-group prediction accuracy: best formula vs best ML "
                   "(each group scored independently; ML on held-out organisms)", fontsize=e_cfg["title_fontsize"],
                   fontweight="bold", pad=e_cfg["title_pad"])

    # Two-part legend: method colors (bars) + category colors (x-tick label
    # text). The category legend previously sat inside a triple-quoted
    # string, which means it was never actually called -- fixed here, and
    # now fully controlled by the panel_e.show_category_legend /
    # category_legend_loc / category_legend_bbox values in CFG.
    method_legend = ax_e.legend(loc=e_cfg["method_legend_loc"], bbox_to_anchor=e_cfg["method_legend_bbox"],
                                fontsize=e_cfg["method_legend_fontsize"], frameon=False)
    ax_e.add_artist(method_legend)

    if e_cfg["show_category_legend"]:
        category_handles = [
            plt.Line2D([], [], color=CATEGORY_COLORS[cat], marker="s", linestyle="",
                       markersize=8, label=cat)
            for cat in ["Eukaryote", "Prokaryote", "Organelle", "Non-cellular"]
        ]
        ax_e.legend(
            handles=category_handles,
            loc=e_cfg["category_legend_loc"],
            bbox_to_anchor=e_cfg["category_legend_bbox"],
            fontsize=e_cfg["category_legend_fontsize"],
            frameon=False,
            title="Group label color",
            title_fontsize=e_cfg["category_legend_fontsize"],
            ncol=e_cfg["category_legend_ncol"],
        )

    style.add_panel_label(ax_e, "E", x=e_cfg["panel_label_x"], y=e_cfg["panel_label_y"])

    # ---- Panel F: summary heatmap, winner Pearson r per base, formula row vs ML row ----
    # (pooled, for a single headline number -- matches how best_*_per_base was chosen)
    # Uses its own nested gridspec so the heatmap itself gets most of the
    # column's width and the colorbar stays a slim strip beside it,
    # rather than matplotlib's default colorbar sizing dominating a
    # heatmap that's only 2 rows x 4 columns of data.
    gs_f = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=gs_bottom[0], width_ratios=f_cfg["width_ratios"], wspace=f_cfg["wspace"]
    )
    ax_f = fig.add_subplot(gs_f[0])
    ax_f_cbar = fig.add_subplot(gs_f[1])

    summary = np.array([
        [_safe_pearson(
            formula_long[(formula_long["Method"] == best_formula_per_base[b]) & (formula_long["Base"] == b)]["Observed"],
            formula_long[(formula_long["Method"] == best_formula_per_base[b]) & (formula_long["Base"] == b)]["Predicted"]
        ) for b in BASES],
        [_safe_pearson(
            ml_long_test[(ml_long_test["Method"] == best_ml_per_base[b]) & (ml_long_test["Base"] == b)]["Observed"],
            ml_long_test[(ml_long_test["Method"] == best_ml_per_base[b]) & (ml_long_test["Base"] == b)]["Predicted"]
        ) for b in BASES],
    ])
    im = ax_f.imshow(summary, cmap=f_cfg["cmap"], vmin=f_cfg["vmin"], vmax=f_cfg["vmax"], aspect="auto")
    ax_f.set_xticks(range(len(BASES)))
    ax_f.set_xticklabels(BASES, fontsize=f_cfg["tick_fontsize"])
    ax_f.set_yticks([0, 1])
    ax_f.set_yticklabels(["Formula", "ML"], fontsize=f_cfg["tick_fontsize"])
    for i in range(2):
        for j in range(len(BASES)):
            ax_f.text(j, i, f"{summary[i,j]:.2f}", ha="center", va="center", fontsize=f_cfg["cell_fontsize"])
    fig.colorbar(im, cax=ax_f_cbar, label="Pearson r")
    ax_f_cbar.tick_params(labelsize=f_cfg["cbar_tick_fontsize"])
    ax_f_cbar.set_ylabel("Pearson r", fontsize=f_cfg["cbar_label_fontsize"])
    ax_f.set_title("Winner summary\n(pooled)", fontsize=f_cfg["title_fontsize"], fontweight="bold")
    style.add_panel_label(ax_f, "F", x=f_cfg["panel_label_x"], y=f_cfg["panel_label_y"])

    # ---- Text panel: which method won per base, placed to F's right ----
    ax_txt = fig.add_subplot(gs_bottom[1])
    ax_txt.axis("off")
    lines = ["Winning method per base (pooled):", ""]
    for b in BASES:
        lines.append(f"  {b}:   Formula = {best_formula_per_base[b]}")
        lines.append(f"       ML = {best_ml_per_base[b]}")
        lines.append("")
    ax_txt.text(0.0, txt_cfg["y_start"], "\n".join(lines), fontsize=txt_cfg["fontsize"], va="top", ha="left",
                transform=ax_txt.transAxes, family="monospace", linespacing=txt_cfg["linespacing"])

    # ---- Panel G: combined CD diagram, formula + ML methods on one rank axis ----
    # Built from the PER-GROUP score tables (the trusted, non-pooled
    # numbers) -- this is the figure's statistical backing for "is the
    # formula-vs-ML difference real, or could it be noise from a few
    # favorable groups?"
    ax_g = fig.add_subplot(gs_bottom[2])
    combined_scores = stats_mod.combine_formula_and_ml_scores(
        formula_scores_per_group, ml_scores_per_group
    )
    c_stat, c_p, c_ranks, c_used, c_total = stats_mod.run_friedman_and_ranks(combined_scores)
    if c_ranks is not None:
        stats_mod.make_combined_cd_diagram(
            c_ranks, c_p, c_used, formula_methods, ml_methods,
            "Combined ranking: formula vs ML\n(all groups & bases)",
            ax=ax_g, compact=True
        )
    else:
        ax_g.axis("off")
        ax_g.text(0.5, 0.5, "Combined CD diagram unavailable\n"
                            "(not enough complete (group, base) blocks\n"
                            "with a score for every method)",
                  ha="center", va="center", fontsize=8.5, color="gray",
                  transform=ax_g.transAxes)

    # Pull whatever legend make_combined_cd_diagram() drew (if any) off
    # the Legend artist itself -- this works whether it was built from
    # real plotted lines or from proxy handles -- then redraw it at the
    # corner set in CFG. 07_stats_cd_diagram.py is never touched.
    old_g_legend = ax_g.get_legend()
    if old_g_legend is not None:
        g_handles = getattr(old_g_legend, "legend_handles", None) or getattr(old_g_legend, "legendHandles", [])
        g_labels = [t.get_text() for t in old_g_legend.get_texts()]
        old_g_legend.remove()
        ax_g.legend(g_handles, g_labels, loc=g_cfg["legend_loc"], bbox_to_anchor=g_cfg["legend_bbox"],
                    fontsize=g_cfg["legend_fontsize"], frameon=False)

    style.add_panel_label(ax_g, "G", x=g_cfg["panel_label_x"], y=g_cfg["panel_label_y"])

    fig.suptitle("Figure: Chromosome nucleotide composition prediction \u2014 "
                 "best formula vs best ML method", fontsize=fig_cfg["suptitle_fontsize"], y=fig_cfg["suptitle_y"])

    style.save_fig(fig, out_path)
    return fig


if __name__ == "__main__":
    paths = importlib.import_module("paths")

    formula_long = pd.read_csv(paths.table_path("formula_results_long.csv"))
    ml_long = pd.read_csv(paths.table_path("ml_results_long.csv"))
    ml_long_test = ml_long[ml_long["Split"] == "test"]

    formula_scores_per_group = pd.read_csv(paths.table_path("formula_scores_per_group.csv"))
    ml_scores_per_group = pd.read_csv(paths.table_path("ml_scores_per_group.csv"))
    formula_scores_pooled = pd.read_csv(paths.table_path("formula_scores_pooled.csv"))
    ml_scores_pooled = pd.read_csv(paths.table_path("ml_scores_pooled.csv"))

    FORMULA_METHODS = ["M1_PlusOnly", "M2_CompOnly", "M3_Average", "M4_Dominant", "M5_Weighted"]
    ML_METHODS = ["RandomForest", "GradientBoost", "SVR", "NeuralNet"]
    # winner per base = lowest MEAN RANK across groups (not the pooled r, which the biggest groups dominate)
    best_formula_per_base = stats_mod.pick_winners_by_mean_rank(formula_scores_per_group, FORMULA_METHODS)
    best_ml_per_base = stats_mod.pick_winners_by_mean_rank(ml_scores_per_group, ML_METHODS)

    make_main_figure(formula_long, ml_long_test, formula_scores_per_group, ml_scores_per_group,
                     best_formula_per_base, best_ml_per_base,
                     paths.figure_path("MainFig_formula_vs_ML"))