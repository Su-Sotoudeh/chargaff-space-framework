"""
06_main_figure.py
-------------------
MAIN FIGURE: Best formula method vs best ML method, across all 4 bases
and all 15 genome groups.

The "best" method per base is chosen by mean rank across groups (see 07_stats_cd_diagram.pick_winners_by_mean_rank)
(one number per method, easy to state in text: "M3_Average was the
best-performing formula method for base A"). But the heatmaps (A, B) and
bar chart (H) are built from the PER-GROUP scores table -- i.e. each
group's own fit quality is never contaminated by another group's data.
The scatter panels (D-G) use predictions generated independently per
group (formula) or trained independently per group (ML). Pooling only
happens for the purpose of choosing one method name and drawing one
combined scatter plot; it never feeds back into how any individual
prediction was produced.

Layout (8 panels), READ TOP-TO-BOTTOM, LEFT-TO-RIGHT (this is the new
order -- the heatmap/CD row used to be at the bottom):

  A: Group2 x formula-method heatmap. Each cell = that method's Pearson r
     averaged over the 4 bases, from the trusted PER-GROUP score table
     (never pooled across groups). Shows the group-name row labels.
  B: Same idea, for the ML methods. Uses the SAME color scale as A (same
     cmap/vmin/vmax), so it carries the one shared colorbar for both A
     and B -- A has none, since it would be redundant. B also has no row
     labels, since A and B are reindexed to the identical group list/
     order, so B's rows already line up with A's.
  C: Combined critical-difference diagram (built in 07_stats_cd_diagram.py).
  D-G: Predicted vs Observed scatter, one per base, in base order A,T,G,C.
     NOTE: panel LETTER and BASE letter are two different things here --
     Panel D shows base A, Panel E shows base T, Panel F shows base G,
     and Panel G shows base C. Easy to misread at a glance since the
     alphabet and the nucleotide letters collide; worth double-checking
     figure captions/text against this mapping.
  H: Grouped bar chart -- per-Group2 Pearson r for the winning formula vs
     winning ML method (winner picked per base, then averaged over the
     4 bases), from the same trusted PER-GROUP table.

===============================================================================
HOW TO ADJUST THIS FIGURE
===============================================================================
Every position, font size, legend location, and spacing value lives in the
CFG dictionary directly below, organized by CFG section rather than by
panel letter (letters keep changing as the layout evolves; section names
don't). Map:

    CFG section        variable name(s)         panel letter(s)
    --------------      -----------------        ---------------
    panel_heatmaps      heat_cfg / ax_a, ax_b     A, B
    panel_cd            cd_cfg   / ax_c           C
    panels_scatter      scatter_cfg / ax (loop)   D, E, F, G
    panel_bar           bar_cfg  / ax_h           H
    panel_positions      (manual overrides)        A, B, C only, so far

Find the section for the panel you want to change, change the value,
re-run. You should never need to edit the plotting code beneath CFG to
move a legend, resize something, or reposition a panel.

NOTE on Panel C's legend: the critical-difference diagram in Panel C is
drawn by stats_mod.make_combined_cd_diagram(), which lives in
07_stats_cd_diagram.py -- NOT in this file. panel_cd only controls this
panel's letter-label position and legend placement (the legend itself is
pulled back out after that function draws it). To change anything else
about the CD diagram's own drawing, that function needs to be shared and
edited directly.

NOTE on panel_positions: setting any of the four values there to a
(x0, y0, width, height) tuple (figure-fraction coordinates) completely
overrides that one panel's automatic position from the GridSpec ratios
above -- use this for final by-hand nudging once the automatic layout is
close but not quite right.
===============================================================================
"""
group_colors = {
    'Chordata': "#1f77b4", 'Viride': "#d62728", 'Insect': "#9467bd",
    'Fungi': "#8c564b", 'Bacteria': "#FFD700", 'Archia': "#7f7f7f",
    'Chloro_plast': "#4682B4", 'Mito': "#0033CC", 'Plasmides': "#8B0000",
    'Virus_DNA+/-': "#4B0082", 'Virus_DNA+': "#FF1493", 'Virus_DNA-': "#A0522D",
    'Virus_RNA+/-': "#800080", 'Virus_RNA+': "#FF8C00", 'Virus_RNA-': "#40E0D0"
}



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

# Biological category color-coding for Panel H's x-axis group labels.
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
# CFG -- every layout / legend / font knob for the figure, grouped by
# CFG-section (see the map in the module docstring above for which
# section controls which panel letter).
# Legend positioning uses matplotlib's normal system:
#   loc          = anchor point on the LEGEND BOX itself (e.g. "lower left"
#                  means the legend's lower-left corner sits at bbox_to_anchor)
#   bbox_to_anchor = (x, y) in AXES-FRACTION coordinates for that panel.
#                  (0,0) = bottom-left of the panel, (1,1) = top-right.
#                  y > 1.0 pushes the legend ABOVE the panel; y < 0 pushes it
#                  below.
# =============================================================================
CFG = {
    "figure": {
        "figsize": (13, 14.5),
        # Row heights, TOP TO BOTTOM as the figure now reads:
        #   row 0 = Panels A, B, C  (heatmaps + CD diagram)
        #   row 1 = Panels D-G      (scatter, one per base)
        #   row 2 = Panel H         (grouped bar chart)
        "height_ratios": [1.65, 1.1, 1.1],
        "hspace": 0.6,
        "wspace": 0.4,
        "top": 0.93,      # top edge of the whole grid, as a fraction of figure height.
        "top_width_ratios": [1.2, 1.3, 1.8],   # A : B : C column widths in the top row
        "top_wspace": 0.3,
        "suptitle_fontsize": 12,
        "suptitle_y": 1.01,
    },

    "panels_scatter": {   # Panels D, E, F, G (bases A, T, G, C respectively -- see docstring)
        "point_size": 9,
        "alpha": 0.55,
        "annotation_fontsize": 6.8,
        "legend_loc": "lower right",
        "legend_bbox": (1.25, 0),
        "legend_fontsize": 6.5,
        "title_fontsize": 10,
        "panel_label_x": -0.25,
        "panel_label_y": 1.18,
    },

    "panel_bar": {   # Panel H
        "bar_width": 0.35,
        "value_label_fontsize": 6.2,
        "value_label_offset": 0.02,
        "xtick_fontsize": 7.5,
        "ylim": (0.2, 1.15),
        "title_fontsize": 10,
        "title_pad": 38,
        "panel_label_x": -0.05,
        "panel_label_y": 1.28,

        "method_legend_loc": "upper left",
        "method_legend_bbox": (0, 1.1),
        "method_legend_fontsize": 7.5,

        "show_category_legend": False,
        "category_legend_loc": "lower left",
        "category_legend_bbox": (0.0, 1.16),
        "category_legend_fontsize": 7.5,
        "category_legend_ncol": 4,
    },

    "panel_heatmaps": {   # Panels A (formula) and B (ML) -- shared cosmetics
        "width_ratios": [4.2, 0.2],   # heatmap : colorbar (used for whichever panel has one)
        "wspace": 0.15,
        "cmap": "RdYlGn",
        "vmin": -1,
        "vmax": 1,
        "cell_fontsize": 6,
        "xtick_fontsize": 7.5,
        "ytick_fontsize": 7,
        "cbar_label_fontsize": 7.5,
        "cbar_tick_fontsize": 7,
        "title_fontsize": 10,
        "formula_label_x": -0.35,   # Panel A's letter position (independent of B's now)
        "formula_label_y": 1.05,
        "ml_label_x": -0.12,        # Panel B's letter position -- smaller magnitude than A's
        "ml_label_y": 1.05,         # since B is a much narrower panel after your manual override
        "formula_title": "Formula methods\n(avg over A,T,G,C)",
        "ml_title": "ML methods\n(avg over A,T,G,C)",
        # A and B are drawn with the identical cmap/vmin/vmax above, which is
        # what makes sharing a single colorbar (B's) valid for both.
        "formula_heatmap_has_colorbar": False,   # Panel A: no colorbar of its own
        "ml_heatmap_has_row_labels": False,      # Panel B: no group-name row labels
    },

    "panel_cd": {   # Panel C
        "panel_label_x": -0.12,
        "panel_label_y": 1.15,
        # Repositions the CD diagram's legend (built inside
        # stats_mod.make_combined_cd_diagram() in 07_stats_cd_diagram.py,
        # not edited here) by pulling it back out and redrawing it here.
        "legend_loc": "upper left",
        "legend_bbox": (0, 0.85),
        "legend_fontsize": 7,
    },

    "panel_positions": {
        # Optional manual overrides for the top row (A, B, C) only. Each
        # value is either None (use the automatic layout from
        # height_ratios / top_width_ratios / top_wspace above) or a
        # 4-tuple (x0, y0, width, height) in FIGURE-fraction coordinates
        # (0,0 = bottom-left corner of the whole figure, 1,1 = top-right),
        # applied with ax.set_position(...) AFTER the automatic layout --
        # this overrides it completely for that one panel, letting you
        # nudge size/position/spacing entirely by hand once the automatic
        # layout is close but not quite right.
        # Example: "formula_heatmap": (0.05, 0.68, 0.26, 0.27)
        "formula_heatmap": (0.16, 0.68, 0.15, 0.25),
        "formula_heatmap_colorbar": None,   # only used if formula_heatmap_has_colorbar is set True above
        "ml_heatmap": (0.33, 0.68, 0.15, 0.25),
        "ml_heatmap_colorbar": (0.49, 0.68, 0.01, 0.25),
        "cd_diagram": (0.58, 0.68, 0.32, 0.25),
    },
}


def _safe_pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    r, _ = pearsonr(x, y)
    return r


def _group_by_method_avg_matrix(scores_per_group, method_order, groups_order, value_col="Pearson_r"):
    """
    Average value_col over whatever Base rows exist (A,T,G,C) for every
    (Group2, Method) pair in a PER-GROUP score table -- i.e. each cell is
    "this method's mean Pearson r for this group, across the 4 bases",
    computed only from that group's own scores (never pooled across
    groups). Returns a Group2 x Method DataFrame ready for imshow.
    pandas' groupby().mean() skips NaN, so a base with an undefined
    Pearson r for that group/method still lets the other 3 bases average
    normally rather than forcing the whole cell to NaN.
    """
    avg = scores_per_group.groupby(["Group2", "Method"])[value_col].mean().reset_index()
    pivot = avg.pivot(index="Group2", columns="Method", values=value_col)
    return pivot.reindex(index=groups_order, columns=method_order)


def _draw_group_method_heatmap(ax, pivot, cfg, cmap, vmin, vmax, ax_cbar=None, show_row_labels=True):
    """
    Draws a Group2 (rows) x Method (columns, labels along the bottom)
    heatmap into ax. If ax_cbar is given, draws a colorbar into it;
    otherwise this panel has no colorbar of its own (relying on a
    neighboring panel's, since both must share the same cmap/vmin/vmax
    for that to be valid). Same NaN convention as
    05_supplementary_figure.py: undefined Pearson r is shown as a gray
    "n/a" cell instead of a silent blank.
    """
    im = ax.imshow(pivot.values, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right", fontsize=cfg["xtick_fontsize"])
    ax.set_yticks(range(pivot.shape[0]))
    if show_row_labels:
        ax.set_yticklabels(pivot.index, fontsize=cfg["ytick_fontsize"])
    else:
        ax.set_yticklabels([])
    for r in range(pivot.shape[0]):
        for c in range(pivot.shape[1]):
            val = pivot.values[r, c]
            if not np.isnan(val):
                ax.text(c, r, f"{val:.2f}", ha="center", va="center", fontsize=cfg["cell_fontsize"])
            else:
                ax.add_patch(plt.Rectangle((c - 0.5, r - 0.5), 1, 1,
                                           facecolor="lightgray", edgecolor="white"))
                ax.text(c, r, "n/a", ha="center", va="center", fontsize=cfg["cell_fontsize"] - 0.5,
                        color="gray", style="italic")
    if ax_cbar is not None:
        fig = ax.figure
        fig.colorbar(im, cax=ax_cbar, label="Pearson r")
        ax_cbar.tick_params(labelsize=cfg["cbar_tick_fontsize"])
        ax_cbar.set_ylabel("Pearson r", fontsize=cfg["cbar_label_fontsize"])
    return im


def _apply_manual_position(ax, manual_pos):
    """
    If manual_pos is a (x0, y0, width, height) tuple in figure-fraction
    coordinates, overrides ax's automatic GridSpec position with it.
    Does nothing if manual_pos is None.
    """
    if manual_pos is not None:
        ax.set_position(manual_pos)


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
    formula_methods / ml_methods : lists of method names for Panels A/B/C
                                  (defaults to the standard 5 formula + 4
                                  ML method names if not provided)
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
    scatter_cfg = cfg["panels_scatter"]
    bar_cfg = cfg["panel_bar"]
    heat_cfg = cfg["panel_heatmaps"]
    cd_cfg = cfg["panel_cd"]
    pos_cfg = cfg["panel_positions"]

    fig = plt.figure(figsize=fig_cfg["figsize"])
    gs = gridspec.GridSpec(3, 4, height_ratios=fig_cfg["height_ratios"],
                           hspace=fig_cfg["hspace"], wspace=fig_cfg["wspace"], top=fig_cfg["top"])

    # ============================================================
    # TOP ROW: Panel A (formula heatmap), Panel B (ML heatmap), Panel C (CD diagram)
    # ============================================================
    gs_top = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=gs[0, :],
        width_ratios=fig_cfg["top_width_ratios"], wspace=fig_cfg["top_wspace"]
    )

    # Both heatmaps are reindexed to this SAME group list/order (rather than
    # each computing its own independently) so Panel B's rows line up
    # exactly with Panel A's -- this is what makes it safe to drop B's row
    # labels below. A group missing from one method's table becomes an
    # "n/a" row there, same convention as 05_supplementary_figure.py.
    groups_shared = [g for g in GROUP_ORDER
                     if g in set(formula_scores_per_group["Group2"].unique())
                     | set(ml_scores_per_group["Group2"].unique())]

    # ---- Panel A: Group x Formula-method heatmap ----
    if heat_cfg["formula_heatmap_has_colorbar"]:
        gs_a = gridspec.GridSpecFromSubplotSpec(
            1, 2, subplot_spec=gs_top[0], width_ratios=heat_cfg["width_ratios"], wspace=heat_cfg["wspace"]
        )
        ax_a = fig.add_subplot(gs_a[0])
        ax_a_cbar = fig.add_subplot(gs_a[1])
    else:
        ax_a = fig.add_subplot(gs_top[0])
        ax_a_cbar = None

    formula_matrix = _group_by_method_avg_matrix(formula_scores_per_group, formula_methods, groups_shared)
    _draw_group_method_heatmap(ax_a, formula_matrix, heat_cfg, heat_cfg["cmap"], heat_cfg["vmin"], heat_cfg["vmax"],
                               ax_cbar=ax_a_cbar, show_row_labels=True)
    ax_a.set_title(heat_cfg["formula_title"], fontsize=heat_cfg["title_fontsize"], fontweight="bold")
    style.add_panel_label(ax_a, "A", x=heat_cfg["formula_label_x"], y=heat_cfg["formula_label_y"])
    _apply_manual_position(ax_a, pos_cfg["formula_heatmap"])
    if ax_a_cbar is not None:
        _apply_manual_position(ax_a_cbar, pos_cfg.get("formula_heatmap_colorbar"))

    # ---- Panel B: Group x ML-method heatmap (carries the shared colorbar) ----
    gs_b = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=gs_top[1], width_ratios=heat_cfg["width_ratios"], wspace=heat_cfg["wspace"]
    )
    ax_b = fig.add_subplot(gs_b[0])
    ax_b_cbar = fig.add_subplot(gs_b[1])
    ml_matrix = _group_by_method_avg_matrix(ml_scores_per_group, ml_methods, groups_shared)
    _draw_group_method_heatmap(ax_b, ml_matrix, heat_cfg, heat_cfg["cmap"], heat_cfg["vmin"], heat_cfg["vmax"],
                               ax_cbar=ax_b_cbar, show_row_labels=heat_cfg["ml_heatmap_has_row_labels"])
    ax_b.set_title(heat_cfg["ml_title"], fontsize=heat_cfg["title_fontsize"], fontweight="bold")
    style.add_panel_label(ax_b, "B", x=heat_cfg["ml_label_x"], y=heat_cfg["ml_label_y"])
    _apply_manual_position(ax_b, pos_cfg["ml_heatmap"])
    _apply_manual_position(ax_b_cbar, pos_cfg["ml_heatmap_colorbar"])

    # ---- Panel C: combined CD diagram, formula + ML methods on one rank axis ----
    # Built from the PER-GROUP score tables (the trusted, non-pooled
    # numbers) -- this is the figure's statistical backing for "is the
    # formula-vs-ML difference real, or could it be noise from a few
    # favorable groups?"
    ax_c = fig.add_subplot(gs_top[2])
    combined_scores = stats_mod.combine_formula_and_ml_scores(
        formula_scores_per_group, ml_scores_per_group
    )
    c_stat, c_p, c_ranks, c_used, c_total = stats_mod.run_friedman_and_ranks(combined_scores)
    if c_ranks is not None:
        stats_mod.make_combined_cd_diagram(
            c_ranks, c_p, c_used, formula_methods, ml_methods,
            "Combined ranking: formula vs ML\n(all groups & bases)",
            ax=ax_c, compact=True
        )
    else:
        ax_c.axis("off")
        ax_c.text(0.5, 0.5, "Combined CD diagram unavailable\n"
                            "(not enough complete (group, base) blocks\n"
                            "with a score for every method)",
                  ha="center", va="center", fontsize=8.5, color="gray",
                  transform=ax_c.transAxes)

    # Pull whatever legend make_combined_cd_diagram() drew (if any) off
    # the Legend artist itself -- works whether it was built from real
    # plotted lines or from proxy handles -- then redraw it at the corner
    # set in CFG. 07_stats_cd_diagram.py is never touched.
    old_c_legend = ax_c.get_legend()
    if old_c_legend is not None:
        c_handles = getattr(old_c_legend, "legend_handles", None) or getattr(old_c_legend, "legendHandles", [])
        c_labels = [t.get_text() for t in old_c_legend.get_texts()]
        old_c_legend.remove()
        ax_c.legend(c_handles, c_labels, loc=cd_cfg["legend_loc"], bbox_to_anchor=cd_cfg["legend_bbox"],
                    fontsize=cd_cfg["legend_fontsize"], frameon=False)

    style.add_panel_label(ax_c, "C", x=cd_cfg["panel_label_x"], y=cd_cfg["panel_label_y"])
    _apply_manual_position(ax_c, pos_cfg["cd_diagram"])

    # ============================================================
    # MIDDLE ROW: Panels D, E, F, G -- scatter per base (pooled FOR DISPLAY ONLY)
    # ============================================================
    panel_labels = iter(["D", "E", "F", "G"])
    for j, base in enumerate(BASES):
        ax = fig.add_subplot(gs[1, j])

        f_method = best_formula_per_base[base]
        m_method = best_ml_per_base[base]

        f_sub = formula_long[(formula_long["Method"] == f_method) & (formula_long["Base"] == base)]
        m_sub = ml_long_test[(ml_long_test["Method"] == m_method) & (ml_long_test["Base"] == base)]

        ax.scatter(f_sub["Observed"], f_sub["Predicted"], s=scatter_cfg["point_size"], alpha=scatter_cfg["alpha"],
                   color=style.METHOD_TYPE_COLORS["formula"], label=f"Formula: {f_method}",
                   linewidths=0, rasterized=True)
        ax.scatter(m_sub["Observed"], m_sub["Predicted"], s=scatter_cfg["point_size"], alpha=scatter_cfg["alpha"],
                   color=style.METHOD_TYPE_COLORS["ml"], label=f"ML: {m_method}",
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
                transform=ax.transAxes, fontsize=scatter_cfg["annotation_fontsize"], va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, edgecolor="lightgray"))

        ax.set_xlabel(f"Observed {base} (%)")
        ax.set_ylabel(f"Predicted {base} (%)")
        ax.set_title(f"Base {base}", fontsize=scatter_cfg["title_fontsize"], fontweight="bold")
        if j == 0:
            ax.legend(loc=scatter_cfg["legend_loc"], bbox_to_anchor=scatter_cfg["legend_bbox"],
                      fontsize=scatter_cfg["legend_fontsize"], frameon=False)
        style.add_panel_label(ax, next(panel_labels), x=scatter_cfg["panel_label_x"], y=scatter_cfg["panel_label_y"])

    # ============================================================
    # BOTTOM ROW: Panel H -- grouped bars, PER-GROUP Pearson r (trusted table, not pooled)
    # ============================================================
    ax_h = fig.add_subplot(gs[2, :])

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
    width = bar_cfg["bar_width"]

    formula_vals = [f_avg.get(g, np.nan) for g in groups_present]
    ml_vals = [m_avg.get(g, np.nan) for g in groups_present]

    bars_f = ax_h.bar(x - width/2, formula_vals, width, label="Best formula (avg over A,T,G,C)",
                      color=style.METHOD_TYPE_COLORS["formula"])
    bars_m = ax_h.bar(x + width/2, ml_vals, width, label="Best ML (avg over A,T,G,C)",
                      color=style.METHOD_TYPE_COLORS["ml"])
    ax_h.axhline(0, color="black", linewidth=0.6)

    for bar, val in zip(bars_f, formula_vals):
        if np.isnan(val):
            continue
        va = "bottom" if val >= 0 else "top"
        offset = bar_cfg["value_label_offset"] if val >= 0 else -bar_cfg["value_label_offset"]
        ax_h.text(bar.get_x() + bar.get_width() / 2, val + offset, f"{val:.2f}",
                  rotation=90, ha="center", va=va, fontsize=bar_cfg["value_label_fontsize"])
    for bar, val in zip(bars_m, ml_vals):
        if np.isnan(val):
            continue
        va = "bottom" if val >= 0 else "top"
        offset = bar_cfg["value_label_offset"] if val >= 0 else -bar_cfg["value_label_offset"]
        ax_h.text(bar.get_x() + bar.get_width() / 2, val + offset, f"{val:.2f}",
                  rotation=90, ha="center", va=va, fontsize=bar_cfg["value_label_fontsize"])

    ax_h.set_xticks(x)
    ax_h.set_xticklabels(groups_present, rotation=45, ha="right", fontsize=bar_cfg["xtick_fontsize"])
    for tick_label, g in zip(ax_h.get_xticklabels(), groups_present):
        cat = GROUP_CATEGORY.get(g, "Eukaryote")
        tick_label.set_color(CATEGORY_COLORS[cat])
        tick_label.set_fontweight("bold")

    ax_h.set_ylabel("Mean Pearson r (across A,T,G,C)")
    ax_h.set_ylim(bar_cfg["ylim"])
    ax_h.set_title("Per-group prediction accuracy: best formula vs best ML "
                   "(each group scored independently; ML on held-out organisms)", fontsize=bar_cfg["title_fontsize"],
                   fontweight="bold", pad=bar_cfg["title_pad"])

    method_legend = ax_h.legend(loc=bar_cfg["method_legend_loc"], bbox_to_anchor=bar_cfg["method_legend_bbox"],
                                fontsize=bar_cfg["method_legend_fontsize"], frameon=False)
    ax_h.add_artist(method_legend)

    if bar_cfg["show_category_legend"]:
        category_handles = [
            plt.Line2D([], [], color=CATEGORY_COLORS[cat], marker="s", linestyle="",
                       markersize=8, label=cat)
            for cat in ["Eukaryote", "Prokaryote", "Organelle", "Non-cellular"]
        ]
        ax_h.legend(
            handles=category_handles,
            loc=bar_cfg["category_legend_loc"],
            bbox_to_anchor=bar_cfg["category_legend_bbox"],
            fontsize=bar_cfg["category_legend_fontsize"],
            frameon=False,
            title="Group label color",
            title_fontsize=bar_cfg["category_legend_fontsize"],
            ncol=bar_cfg["category_legend_ncol"],
        )

    style.add_panel_label(ax_h, "H", x=bar_cfg["panel_label_x"], y=bar_cfg["panel_label_y"])

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