"""
Chargaff Fig 6 -- full pipeline, FROM ZERO every run.

* No existing PNG/SVG/TIFF is reused: old outputs are wiped at start and every
  panel is rendered fresh.
* The final SVG is assembled directly (one real SVG document):
    - bar panel E, panel letters, F-strip headers/swatches/text = true vector
    - 3D panels A-D and the 16 F cells = high-res PNG <image> (Plotly 3D is
      WebGL, so it can only be rasterised; rendered at RENDER_SCALE_3D x)
* Raster composite (PNG + 600-dpi TIFF) is also built, from the same fresh panels.

Outputs
-------
fig6_panels/                 individual fresh PNG panels (3D at RENDER_SCALE_3D x)
publication_exports/
    chargaff_full_figure.svg        pure-SVG assembled figure
    chargaff_full_figure.png/.tiff  raster composite (EXPORT_SCALE x)
    chargaff_full_figure_outlined.svg / .pdf   (if Inkscape is installed)
    panel_E_bar.svg                 true-vector bar panel on its own
    fig6_group_stats.csv
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)
import base64
import io
import os
import re
import shutil
import subprocess
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from PIL import Image, ImageDraw, ImageFont

Image.MAX_IMAGE_PIXELS = None

# master table: Dataframes/<master>.csv (override with CSPR_MASTER_CSV or CHARGAFF_CSV, see config.py)
file_path = str(config.require(config.MASTER_CSV, 'Copy the master table into Dataframes/ (see README.md).'))

# DejaVu fonts ship with matplotlib, so they exist on every OS (PIL cannot find them by bare name on Windows/macOS)
try:
    import matplotlib
    _FONT_DIR = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
except Exception:
    _FONT_DIR = ""

# ============================================================================
# CONFIGURATION
# ============================================================================

USE_SOFTWARE_RENDERING = True
CLEAN_START            = True    # delete old files inside PANEL_DIR and OUT_DIR first

M_MIN = 0.971716
M_MAX = 1.028284
Z_COLUMN = "Z"

# ---- 3D geometry (shared) --------------------------------------------------
GRID_RESOLUTION            = 50
VOLUME_OPACITY             = 0.10
VOLUME_SURFACE_COUNT       = 15
VOLUME_COLOR               = "green"

WIREFRAME_EDGE_COLOR       = "black"
WIREFRAME_EDGE_WIDTH       = 2
WIREFRAME_EDGE_RESOLUTION  = 300

CONTOUR_RING_COLOR         = "rgba(255,255,255,0.65)"
CONTOUR_RING_WIDTH         = 1.5
CONTOUR_RING_RESOLUTION    = 300
CONTOUR_Z_LEVELS_UPPER     = [99, 95, 90, 82, 75, 68, 57]
CONTOUR_Z_LEVELS_LOWER     = [1, 5, 10, 18, 25, 32, 43]

EQUATOR_RADIUS             = 35
EQUATOR_RESOLUTION         = 400
EQUATOR_COLOR_OUTER        = "gold"
EQUATOR_WIDTH_OUTER        = 4
EQUATOR_COLOR_INNER        = "darkgoldenrod"
EQUATOR_WIDTH_INNER        = 1.5
EQUATOR_Z_OFFSET           = 1.2

AXIS_LINE_COLOR            = "black"
AXIS_LINE_WIDTH            = 3
AXIS_LINE_Z_OFFSET         = 0.9

CORNER_LABEL_FONT_SIZE     = 45   # A/T/C/G and +/- labels in the big panels A-D
CORNER_LABEL_FONT_SIZE_STRIP = 18  # same labels in the small F-strip panels
CORNER_LABEL_COLOR         = "black"
CORNER_LABEL_Z_OFFSET      = -0.5
POLE_LABEL_TOP_OFFSET      = 10
POLE_LABEL_BOTTOM_OFFSET   = 10

SCENE_AXIS_PAD             = 5
SCENE_Z_PAD_TOP            = 10
SCENE_Z_PAD_BOTTOM         = 10
CAMERA_EYE                 = dict(x=1.15, y=1.15, z=0.18)

# ---- Big panels A, B, C, D -------------------------------------------------
PANEL_SIZE         = 1400
PANEL_MARGIN       = dict(l=0, r=0, t=0, b=0)
MARKER_SIZE        = 5
MARKER_OPACITY     = 1.0

# ---- Scales ----------------------------------------------------------------
# EXPORT_SCALE    : layout unit of the figure (SVG user units / raster composite
#                   pixels).  Font sizes and offsets below are in these units.
# RENDER_SCALE_3D : pixel density of the 3D PNGs that get embedded in the SVG.
#                   Higher = sharper zoom, bigger file, slower render.
EXPORT_SCALE       = 2
RENDER_SCALE_3D    = 2     # 2 is already ~1250 dpi at 18 cm figure width; 4 only inflates the file

# ---- PNG size reduction for the images embedded in the SVG -----------------
# Lossless recompression is always applied (Pillow optimize; plus oxipng if installed).
# PNG_PALETTE_COLORS: None = fully lossless. 256 (or 128) = adaptive palette, no dithering;
# usually visually identical for flat-colour scatter panels and typically 2-3x smaller,
# but CHECK the translucent green volume for banding before using it for submission.
PNG_PALETTE_COLORS = None

# ---- Small panel render size for the F strip ------------------------------
BAR_PANEL_WIDTH    = 1600
SMALL_PANEL_SIZE   = (2 * PANEL_SIZE + BAR_PANEL_WIDTH) // 8   # = 550
STRIP_MARKER_SIZE  = 2

# ---- Sampling caps ---------------------------------------------------------
SAMPLE_LIMIT_MAIN      = None
SAMPLE_LIMIT_ORG_VIRUS = 5000
SAMPLE_RANDOM_STATE    = 1

# ---- Bar panel E -----------------------------------------------------------
BAR_PANEL_HEIGHT          = 2 * PANEL_SIZE
BAR_CENTER                = 50
BAR_X_AXIS_RANGE          = [-55, 165]
BAR_LABEL_GAP             = 2
BAR_LABEL_FONT_SIZE       = 50
BAR_GROUP_LABEL_FONT_SIZE = 65
BAR_MARGIN                = dict(l=10, r=10, t=20, b=40)

# ---- Panel labels A-F ------------------------------------------------------
PANEL_LABEL_FONT_SIZE = 90
PANEL_LABEL_FONT_PATH = os.path.join(_FONT_DIR, "DejaVuSans-Bold.ttf")
PANEL_LABEL_OFFSET    = (40, 5)

# ---- F strip --------------------------------------------------------------
GROUPS_FIG6 = [
    'Mammals', 'Chordata', 'Plantae', 'Arthropoda', 'Fungi', 'Bacteria', 'Archaea', 'Chloroplast',
    'Plasmids', 'Mitochondria',
    'Virus_DNA+/-', 'Virus_DNA+', 'Virus_DNA-', 'Virus_RNA+/-', 'Virus_RNA+', 'Virus_RNA-',
]
STRIP_COLUMNS            = 8
STRIP_HEADER_HEIGHT      = 200
HEADER_NAME_FONT_SIZE    = 95
HEADER_STATS_FONT_SIZE   = 75
HEADER_SWATCH_SIZE       = 90
HEADER_X_OFFSET          = 30
HEADER_FONT_PATH         = os.path.join(_FONT_DIR, "DejaVuSans-Bold.ttf")
HEADER_STATS_FONT_PATH   = os.path.join(_FONT_DIR, "DejaVuSans.ttf")
GAP_BEFORE_STRIP         = 130

# ---- Output ----------------------------------------------------------------
PANEL_DIR   = str(config.output_dir(6, "fig6_panels"))
OUT_DIR     = str(config.output_dir(6, "publication_exports"))
FULL_NAME   = "chargaff_full_figure"
STATS_CSV   = os.path.join(OUT_DIR, "fig6_group_stats.csv")

EMBED_IMAGES_IN_SVG     = True    # False -> SVG links PNGs in PANEL_DIR (small file)
EXPORT_RASTER_COMPOSITE = True    # PNG + TIFF composite
TIFF_DPI                = 600
INKSCAPE_OUTLINED       = True    # if inkscape is on PATH: outlined-text SVG + PDF

# SVG text font stack (viewer picks the first installed). Outline text for submission.
SVG_FONT_FAMILY = "DejaVu Sans, Verdana, Arial, sans-serif"
SVG_ASCENT      = 0.928           # DejaVu Sans ascent/em (PIL text y = ascender line)

# ---- Group lists ----------------------------------------------------------
groups_ABCD = [
    'Chordata', 'Plantae', 'Arthropoda', 'Fungi', 'Bacteria', 'Archaea',
    'Chloroplast', 'Mitochondria', 'Plasmids',
    'Virus_DNA+/-', 'Virus_DNA+', 'Virus_DNA-',
    'Virus_RNA+/-', 'Virus_RNA+', 'Virus_RNA-',
]
groups_BAR   = ['Mammals'] + groups_ABCD
groups_STATS = list(dict.fromkeys(groups_BAR + GROUPS_FIG6))

group_colors = {
    'Mammals':       "#17becf",
    'Chordata':      "#1f77b4",
    'Plantae':       "#d62728",
    'Arthropoda':    "#9467bd",
    'Fungi':         "#8c564b",
    'Bacteria':      "#FFD700",
    'Archaea':       "#7f7f7f",
    'Chloroplast':   "#4682B4",
    'Mitochondria':  "#0033CC",
    'Plasmids':      "#8B0000",
    'Virus_DNA+/-':  "#4B0082",
    'Virus_DNA+':    "#FF1493",
    'Virus_DNA-':    "#A0522D",
    'Virus_RNA+/-':  "#800080",
    'Virus_RNA+':    "#FF8C00",
    'Virus_RNA-':    "#40E0D0",
}

# ============================================================================
# FIGURE LAYOUT (in EXPORT_SCALE units; shared by the SVG and raster composite)
# ============================================================================

PANEL_PX  = PANEL_SIZE * EXPORT_SCALE
BAR_W_PX  = BAR_PANEL_WIDTH * EXPORT_SCALE
BAR_H_PX  = BAR_PANEL_HEIGHT * EXPORT_SCALE
LEFT_W    = 2 * PANEL_PX
TOTAL_W   = LEFT_W + BAR_W_PX
MAIN_H    = 2 * PANEL_PX

CELL_PX   = TOTAL_W // STRIP_COLUMNS
N_ROWS    = -(-len(GROUPS_FIG6) // STRIP_COLUMNS)
ROW_H     = STRIP_HEADER_HEIGHT + CELL_PX
STRIP_W   = STRIP_COLUMNS * CELL_PX
STRIP_X0  = (TOTAL_W - STRIP_W) // 2
STRIP_Y0  = MAIN_H + GAP_BEFORE_STRIP
TOTAL_H   = STRIP_Y0 + N_ROWS * ROW_H

PANEL_POS = {"A": (0, 0), "B": (PANEL_PX, 0), "C": (0, PANEL_PX), "D": (PANEL_PX, PANEL_PX)}

if SMALL_PANEL_SIZE * EXPORT_SCALE != CELL_PX:
    print(f"WARNING: strip panel render size ({SMALL_PANEL_SIZE * EXPORT_SCALE}) != "
          f"strip cell ({CELL_PX}); cells will be slightly rescaled.")

# ============================================================================
# SETUP
# ============================================================================

if USE_SOFTWARE_RENDERING:
    try:
        pio.kaleido.scope.chromium_args = tuple(pio.kaleido.scope.chromium_args) + (
            "--disable-gpu", "--use-gl=swiftshader",
            "--enable-webgl", "--ignore-gpu-blocklist",
        )
    except AttributeError:
        print("kaleido.scope.chromium_args not available — kaleido>=1.0? "
              "Consider `pip install kaleido==0.2.1`.")


def clean_outputs():
    """Remove old figure files so NOTHING from a previous run can be reused.
    Only touches image/vector files inside PANEL_DIR and OUT_DIR."""
    for d in (PANEL_DIR, OUT_DIR):
        if os.path.isdir(d):
            for fn in os.listdir(d):
                if fn.lower().endswith((".png", ".svg", ".pdf", ".tif", ".tiff", ".csv")):
                    os.remove(os.path.join(d, fn))


if CLEAN_START:
    clean_outputs()
os.makedirs(PANEL_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================================
# LOAD DATA
# ============================================================================

print("\n=========================================================")
print("LOADING DATA")
print("=========================================================")

df = pd.read_csv(file_path, low_memory=False)
print(f"Loaded {len(df):,} rows, {len(df.columns)} columns")

if Z_COLUMN not in df.columns:
    raise ValueError(f"Z column '{Z_COLUMN}' not found in dataframe")


def virus_subgroup(row):
    g = str(row["Group"]).strip()
    if g != "Viruses":
        return g
    gt = str(row["Genome_Type"]).upper()
    cd = str(row["Coding"]).replace(" ", "")
    if "RNA" in gt:
        return "Virus_RNA+" if cd == "+" else "Virus_RNA-" if cd == "-" else "Virus_RNA+/-"
    if "DNA" in gt:
        return "Virus_DNA+" if cd == "+" else "Virus_DNA-" if cd == "-" else "Virus_DNA+/-"
    return "Unknown"


df["Group2"] = df.apply(virus_subgroup, axis=1)

print("\nGroup2 counts:")
print(df["Group2"].value_counts(dropna=False).to_string())

mam_n = int((df["Group2"] == "Mammals").sum())
print(f"\n>>> Mammals rows in Group2: {mam_n:,}")
if mam_n == 0:
    print("    WARNING: no Mammals rows — the bar plot will have no Mammals entry.")


def sample_limit_for(name):
    if name in ('Chloroplast', 'Mitochondria', 'Plasmids', 'Virus_DNA+/-'):
        return SAMPLE_LIMIT_ORG_VIRUS
    return SAMPLE_LIMIT_MAIN


# ============================================================================
# GRID + SCENE
# ============================================================================

X_vals = np.linspace(0, 100, GRID_RESOLUTION)
Xg, Yg, Zg = np.meshgrid(X_vals, X_vals, X_vals, indexing='ij')
dxg = (Xg - 50) / 50
dyg = (Yg - 50) / 50
dzg = (Zg - 50) / 50
Mg  = 1 + dzg * np.sqrt(dxg**2 + dyg**2)

scene_config = dict(
    xaxis=dict(visible=False, showbackground=False, showgrid=False, zeroline=False,
               range=[-SCENE_AXIS_PAD, 100 + SCENE_AXIS_PAD]),
    yaxis=dict(visible=False, showbackground=False, showgrid=False, zeroline=False,
               range=[-SCENE_AXIS_PAD, 100 + SCENE_AXIS_PAD]),
    zaxis=dict(visible=False, showbackground=False, showgrid=False, zeroline=False,
               range=[-SCENE_Z_PAD_BOTTOM, 100 + SCENE_Z_PAD_TOP]),
    aspectmode='cube',
    camera=dict(eye=CAMERA_EYE),
    bgcolor='rgba(0,0,0,0)'
)


# ============================================================================
# SPACE TRACES
# ============================================================================

def build_space_traces(label_font_size=CORNER_LABEL_FONT_SIZE):
    traces = []

    traces.append(go.Volume(
        x=Xg.flatten(), y=Yg.flatten(), z=Zg.flatten(), value=Mg.flatten(),
        isomin=M_MIN, isomax=M_MAX, opacity=VOLUME_OPACITY,
        surface_count=VOLUME_SURFACE_COUNT,
        colorscale=[[0, VOLUME_COLOR], [1, VOLUME_COLOR]],
        showscale=False, showlegend=False, hoverinfo='skip'
    ))

    n_side = WIREFRAME_EDGE_RESOLUTION
    x1 = np.linspace(0, 100, n_side); y1 = np.zeros(n_side)
    x2 = np.full(n_side, 100.0);     y2 = np.linspace(0, 100, n_side)
    x3 = np.linspace(100, 0, n_side); y3 = np.full(n_side, 100.0)
    x4 = np.full(n_side, 0.0);        y4 = np.linspace(100, 0, n_side)
    edge_x = np.concatenate([x1, x2, x3, x4])
    edge_y = np.concatenate([y1, y2, y3, y4])

    dx_edge = (edge_x - 50) / 50
    dy_edge = (edge_y - 50) / 50
    rho_edge = np.maximum(np.sqrt(dx_edge**2 + dy_edge**2), 1e-12)
    delta   = M_MAX - 1.0
    delta_z = 50 * (delta / rho_edge)

    for ez in (50 + delta_z, 50 - delta_z):
        traces.append(go.Scatter3d(
            x=edge_x, y=edge_y, z=ez, mode='lines',
            line=dict(color=WIREFRAME_EDGE_COLOR, width=WIREFRAME_EDGE_WIDTH),
            showlegend=False, hoverinfo='skip'
        ))

    theta = np.linspace(0, 2 * np.pi, CONTOUR_RING_RESOLUTION)

    def surface_radius_from_z(z):
        dz = abs((z - 50) / 50)
        if dz < 1e-12:
            return None
        return 50 * (delta / dz)

    for zc in CONTOUR_Z_LEVELS_UPPER + CONTOUR_Z_LEVELS_LOWER:
        r_mesh = surface_radius_from_z(zc)
        if r_mesh is not None and r_mesh <= 50:
            traces.append(go.Scatter3d(
                x=50 + r_mesh * np.cos(theta),
                y=50 + r_mesh * np.sin(theta),
                z=np.full_like(theta, zc),
                mode='lines',
                line=dict(color=CONTOUR_RING_COLOR, width=CONTOUR_RING_WIDTH),
                showlegend=False, hoverinfo='skip'
            ))

    theta_gold = np.linspace(0, 2 * np.pi, EQUATOR_RESOLUTION)
    circle_x = 50 + EQUATOR_RADIUS * np.cos(theta_gold)
    circle_y = 50 + EQUATOR_RADIUS * np.sin(theta_gold)
    circle_z = np.full_like(theta_gold, 50 + EQUATOR_Z_OFFSET)
    traces.append(go.Scatter3d(
        x=circle_x, y=circle_y, z=circle_z, mode='lines',
        line=dict(color=EQUATOR_COLOR_OUTER, width=EQUATOR_WIDTH_OUTER),
        showlegend=False, hoverinfo='skip'
    ))
    traces.append(go.Scatter3d(
        x=circle_x, y=circle_y, z=circle_z + 0.2, mode='lines',
        line=dict(color=EQUATOR_COLOR_INNER, width=EQUATOR_WIDTH_INNER),
        showlegend=False, hoverinfo='skip'
    ))

    axis_z = 50 + AXIS_LINE_Z_OFFSET
    for ax, ay, az in (
            ([0, 100], [50, 50], [axis_z, axis_z]),
            ([50, 50], [0, 100], [axis_z, axis_z]),
            ([50, 50], [50, 50], [0, 100])
    ):
        traces.append(go.Scatter3d(
            x=ax, y=ay, z=az, mode='lines',
            line=dict(color=AXIS_LINE_COLOR, width=AXIS_LINE_WIDTH),
            showlegend=False, hoverinfo='skip'
        ))

    labels = [
        (100, 50, 50 - CORNER_LABEL_Z_OFFSET, "A"),
        (0,   50, 50 - CORNER_LABEL_Z_OFFSET, "T"),
        (50, 100, 50 - CORNER_LABEL_Z_OFFSET, "C"),
        (50,   0, 50 - CORNER_LABEL_Z_OFFSET, "G"),
        (50,  50, 100 + POLE_LABEL_TOP_OFFSET, "+"),
        (50,  50, 0   - POLE_LABEL_BOTTOM_OFFSET, "−"),
    ]
    for x, y, z, text in labels:
        traces.append(go.Scatter3d(
            x=[x], y=[y], z=[z], mode='text', text=[text],
            textfont=dict(size=label_font_size, color=CORNER_LABEL_COLOR),
            showlegend=False, hoverinfo='skip'
        ))

    return traces


space_traces       = build_space_traces(CORNER_LABEL_FONT_SIZE)        # panels A-D
space_traces_strip = build_space_traces(CORNER_LABEL_FONT_SIZE_STRIP)  # F strip


# ============================================================================
# DATA COMPUTATION
# ============================================================================

def compute_points(dfg):
    if '-' in str(dfg["Group2"].iloc[0]):
        A = pd.to_numeric(dfg["annot_C_A"], errors='coerce')
        T = pd.to_numeric(dfg["annot_C_T"], errors='coerce')
        G = pd.to_numeric(dfg["annot_C_G"], errors='coerce')
        C = pd.to_numeric(dfg["annot_C_C"], errors='coerce')
    else:
        A = pd.to_numeric(dfg["annot_P_A"], errors='coerce')
        T = pd.to_numeric(dfg["annot_P_T"], errors='coerce')
        G = pd.to_numeric(dfg["annot_P_G"], errors='coerce')
        C = pd.to_numeric(dfg["annot_P_C"], errors='coerce')

    P = pd.to_numeric(dfg[Z_COLUMN], errors='coerce')

    X = 100 * A / (A + T)
    Y = 100 * G / (G + C)
    Z = P

    dx = (X - 50) / 50
    dy = (Y - 50) / 50
    dz = (Z - 50) / 50
    M = 1 + dz * np.sqrt(dx**2 + dy**2)

    valid = ~(np.isnan(X) | np.isnan(Y) | np.isnan(Z) | np.isnan(M))
    return X[valid].values, Y[valid].values, Z[valid].values, M[valid].values


data_types_panel_a = {
    'annot_P': {'A': 'annot_P_A', 'T': 'annot_P_T', 'G': 'annot_P_G', 'C': 'annot_P_C', 'Z_fixed': 100},
    'annot_C': {'A': 'annot_C_A', 'T': 'annot_C_T', 'G': 'annot_C_G', 'C': 'annot_C_C', 'Z_fixed': 0},
    'Unanot':  {'A': 'Unanot_A',  'T': 'Unanot_T',  'G': 'Unanot_G',  'C': 'Unanot_C',  'Z_fixed': 50},
}


def compute_points_panel_a(dfg, dtype_info):
    A = pd.to_numeric(dfg[dtype_info['A']], errors='coerce')
    T = pd.to_numeric(dfg[dtype_info['T']], errors='coerce')
    G = pd.to_numeric(dfg[dtype_info['G']], errors='coerce')
    C = pd.to_numeric(dfg[dtype_info['C']], errors='coerce')

    AT = A + T
    GC = G + C
    X = 100 * A / AT
    Y = 100 * G / GC
    Z = np.full(len(dfg), dtype_info['Z_fixed'], dtype=float)

    dx = (X - 50) / 50
    dy = (Y - 50) / 50
    dz = (Z - 50) / 50
    M = 1 + dz * np.sqrt(dx**2 + dy**2)

    valid = np.isfinite(X) & np.isfinite(Y) & np.isfinite(M) & (AT > 0) & (GC > 0)
    return X[valid].values, Y[valid].values, Z[valid], M[valid]


# ============================================================================
# POINTS PER GROUP (computed once, reused by A-D, stats and the F strip)
# ============================================================================

print("\n=========================================================")
print("COMPUTING POINTS AND STATISTICS")
print("=========================================================")

sampled = {}      # group -> sampled dataframe
points  = {}      # group -> (X, Y, Z, M)
stats   = []

for group in groups_STATS:
    dfg = df[df["Group2"] == group]
    limit = sample_limit_for(group)
    if limit is not None and len(dfg) > limit:
        dfg = dfg.sample(n=limit, random_state=SAMPLE_RANDOM_STATE)
    if len(dfg) == 0:
        print(f"  {group:<18} -> 0 rows (skipped)")
        continue

    X, Y, Z, M = compute_points(dfg)
    if len(M) == 0:
        print(f"  {group:<18} -> no valid points (skipped)")
        continue

    sampled[group] = dfg
    points[group]  = (X, Y, Z, M)

    inside = int(np.sum((M >= M_MIN) & (M <= M_MAX)))
    total  = len(M)
    stats.append({
        "Group":          group,
        "Total_samples":  total,
        "Inside_space":   inside,
        "Outside_space":  total - inside,
        "Percent_inside": 100 * inside / total,
    })
    print(f"  {group:<18} -> {total:>10,} pts, {100 * inside / total:6.2f}% inside")

stats_df = pd.DataFrame(stats)
stats_lookup = {row["Group"]: row for _, row in stats_df.iterrows()}

print("\nGroup statistics:")
print(stats_df.sort_values("Percent_inside", ascending=False)
      .to_string(index=False, formatters={"Percent_inside": "{:.2f}".format}))
stats_df.to_csv(STATS_CSV, index=False)


# ============================================================================
# A-D TRACES
# ============================================================================

traces_A, traces_B, traces_C, traces_D = [], [], [], []

for group in groups_ABCD:
    if group not in points:
        continue
    dfg = sampled[group]
    X, Y, Z, M = points[group]

    mask_inside = (M >= M_MIN) & (M <= M_MAX)
    marker_style = dict(size=MARKER_SIZE, color=group_colors[group], opacity=MARKER_OPACITY)

    for _, dtype_info in data_types_panel_a.items():
        required = [dtype_info['A'], dtype_info['T'], dtype_info['G'], dtype_info['C']]
        if not all(c in dfg.columns for c in required):
            continue
        Xa, Ya, Za, _ = compute_points_panel_a(dfg, dtype_info)
        if len(Xa) == 0:
            continue
        traces_A.append(go.Scatter3d(
            x=Xa, y=Ya, z=Za, mode='markers',
            marker=marker_style, showlegend=False, hoverinfo='skip'))

    traces_B.append(go.Scatter3d(
        x=X, y=Y, z=Z, mode='markers',
        marker=marker_style, showlegend=False, hoverinfo='skip'))
    traces_C.append(go.Scatter3d(
        x=X[mask_inside], y=Y[mask_inside], z=Z[mask_inside], mode='markers',
        marker=marker_style, showlegend=False, hoverinfo='skip'))
    traces_D.append(go.Scatter3d(
        x=X[~mask_inside], y=Y[~mask_inside], z=Z[~mask_inside], mode='markers',
        marker=marker_style, showlegend=False, hoverinfo='skip'))


# ============================================================================
# RENDER 3D PANELS (always fresh, PNG at RENDER_SCALE_3D x)
# ============================================================================

def shrink_png(path):
    """Reduce PNG file size. Lossless unless PNG_PALETTE_COLORS is set."""
    before = os.path.getsize(path)
    im = Image.open(path).convert("RGB")
    if PNG_PALETTE_COLORS:
        im = im.quantize(colors=PNG_PALETTE_COLORS, method=Image.Quantize.MEDIANCUT,
                         dither=Image.Dither.NONE)
    im.save(path, format="PNG", optimize=True)
    ox = shutil.which("oxipng")
    if ox:
        try:
            subprocess.run([ox, "-o", "4", "--strip", "all", "-q", path], check=True)
        except Exception as e:
            print(f"  oxipng failed ({e}); kept Pillow-optimized PNG")
    after = os.path.getsize(path)
    print(f"  shrunk {os.path.basename(path)}: {before / 1e6:.2f} -> {after / 1e6:.2f} MB")


def render_3d_png(traces, path, size):
    fig = go.Figure()
    fig.add_traces(space_traces)
    for tr in traces:
        fig.add_trace(tr)
    fig.update_layout(width=size, height=size, scene=scene_config,
                      paper_bgcolor='white', margin=PANEL_MARGIN, showlegend=False)
    fig.write_image(path, format="png", width=size, height=size, scale=RENDER_SCALE_3D)
    print(f"Rendered {path}  ({size * RENDER_SCALE_3D}px)")
    shrink_png(path)
    return path


print("\n=========================================================")
print("RENDERING PANELS A, B, C, D")
print("=========================================================")

PANEL_FILES = {
    "A": render_3d_png(traces_A, os.path.join(PANEL_DIR, "panel_A.png"), PANEL_SIZE),
    "B": render_3d_png(traces_B, os.path.join(PANEL_DIR, "panel_B.png"), PANEL_SIZE),
    "C": render_3d_png(traces_C, os.path.join(PANEL_DIR, "panel_C.png"), PANEL_SIZE),
    "D": render_3d_png(traces_D, os.path.join(PANEL_DIR, "panel_D.png"), PANEL_SIZE),
}


# ============================================================================
# F STRIP PANELS
# ============================================================================

def strip_panel_path(group):
    safe = group.replace("+", "p").replace("-", "m").replace("/", "")
    return os.path.join(PANEL_DIR, f"fig6_{safe}.png")


print("\n=========================================================")
print(f"RENDERING F STRIP PANELS ({len(GROUPS_FIG6)})")
print("=========================================================")

for group in GROUPS_FIG6:
    if group not in points:
        print(f"  {group}: no data, no panel")
        continue
    X, Y, Z, _ = points[group]
    fig = go.Figure()
    fig.add_traces(space_traces_strip)
    fig.add_trace(go.Scatter3d(
        x=X, y=Y, z=Z, mode='markers',
        marker=dict(size=STRIP_MARKER_SIZE, color=group_colors[group], opacity=MARKER_OPACITY),
        showlegend=False, hoverinfo='skip'))
    fig.update_layout(width=SMALL_PANEL_SIZE, height=SMALL_PANEL_SIZE, scene=scene_config,
                      paper_bgcolor='white', margin=PANEL_MARGIN, showlegend=False)
    out = strip_panel_path(group)
    fig.write_image(out, format="png", width=SMALL_PANEL_SIZE, height=SMALL_PANEL_SIZE,
                    scale=RENDER_SCALE_3D)
    print(f"Rendered {out}")
    shrink_png(out)


# ============================================================================
# BAR PANEL E  (2D -> true vector SVG; PNG only for the raster composite)
# ============================================================================

def build_bar_fig(stats_df, width, height):
    bar_df = stats_df[stats_df["Group"].isin(groups_BAR)].copy()
    missing = [g for g in groups_BAR if g not in bar_df["Group"].values]
    if missing:
        print(f"  NOTE: bar plot is missing these groups (no stats): {missing}")
    print(f"  Bar plot groups: {list(bar_df['Group'])}")

    bar_df = bar_df.sort_values("Percent_inside", ascending=False).reset_index(drop=True)
    values = bar_df["Percent_inside"]
    base   = BAR_CENTER - values / 2
    end    = base + values

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=values, y=bar_df["Group"], base=base, orientation='h',
        marker=dict(color=[group_colors[g] for g in bar_df["Group"]],
                    line=dict(color='black', width=1)),
        showlegend=False, hoverinfo='skip'))
    fig.add_trace(go.Scatter(
        x=base - BAR_LABEL_GAP, y=bar_df["Group"], mode='text',
        text=bar_df["Group"], textposition='middle left',
        textfont=dict(size=BAR_GROUP_LABEL_FONT_SIZE, color='black'),
        showlegend=False, hoverinfo='skip'))
    fig.add_trace(go.Scatter(
        x=end + BAR_LABEL_GAP, y=bar_df["Group"], mode='text',
        text=[f"{p:.1f}%  (n={n:,})" for p, n in
              zip(bar_df["Percent_inside"], bar_df["Total_samples"])],
        textposition='middle right',
        textfont=dict(size=BAR_LABEL_FONT_SIZE, color='black'),
        showlegend=False, hoverinfo='skip'))
    fig.update_xaxes(range=BAR_X_AXIS_RANGE, showticklabels=False,
                     showgrid=False, zeroline=False, title_text=None)
    fig.update_yaxes(showticklabels=False, title_text=None)
    fig.update_layout(width=width, height=height,
                      paper_bgcolor='white', plot_bgcolor='white', margin=BAR_MARGIN)
    return fig


print("\n=========================================================")
print("RENDERING BAR PANEL E")
print("=========================================================")

bar_fig = build_bar_fig(stats_df, BAR_PANEL_WIDTH, BAR_PANEL_HEIGHT)
BAR_PNG_FILE = os.path.join(PANEL_DIR, "panel_E_bar.png")
bar_fig.write_image(BAR_PNG_FILE, format="png", width=BAR_PANEL_WIDTH,
                    height=BAR_PANEL_HEIGHT, scale=EXPORT_SCALE)
bar_svg = bar_fig.to_image(format="svg", width=BAR_PANEL_WIDTH,
                           height=BAR_PANEL_HEIGHT, scale=1).decode("utf-8")
with open(os.path.join(OUT_DIR, "panel_E_bar.svg"), "w", encoding="utf-8") as f:
    f.write(bar_svg)
print("Rendered bar panel (PNG + true-vector SVG)")


# ============================================================================
# PURE-SVG ASSEMBLY
# ============================================================================

def _svg_image(path, x, y, w, h, embed, svg_dir):
    if embed:
        with open(path, "rb") as fh:
            href = "data:image/png;base64," + base64.b64encode(fh.read()).decode("ascii")
    else:
        href = os.path.relpath(path, svg_dir).replace(os.sep, "/")
    return (f'<image x="{x}" y="{y}" width="{w}" height="{h}" '
            f'preserveAspectRatio="none" xlink:href="{href}"/>')


def _svg_text(x, y_top, s, size, bold=False, fill="black"):
    """y_top matches PIL's draw.text((x, y_top), ...) -> convert to SVG baseline."""
    return (f'<text x="{x}" y="{y_top + SVG_ASCENT * size:.1f}" '
            f'font-family="{SVG_FONT_FAMILY}" font-size="{size}" '
            f'font-weight="{"bold" if bold else "normal"}" fill="{fill}">{escape(s)}</text>')


def _inline_plotly_svg(svg_text, x, y, target_w, target_h):
    m = re.search(r"<svg\b[^>]*>", svg_text)
    if m is None:
        raise ValueError("bar SVG does not look like an SVG document")
    root = m.group(0)
    w = float(re.search(r'\swidth="([\d.]+)', root).group(1))
    h = float(re.search(r'\sheight="([\d.]+)', root).group(1))
    inner = svg_text[m.end(): svg_text.rindex("</svg>")]
    return (f'<g transform="translate({x},{y}) scale({target_w / w:.6f},{target_h / h:.6f})">'
            f"{inner}</g>")


def assemble_svg(out_path):
    svg_dir = os.path.dirname(os.path.abspath(out_path))
    ox, oy = PANEL_LABEL_OFFSET
    embed = EMBED_IMAGES_IN_SVG

    with open(out_path, "w", encoding="utf-8") as f:
        w = f.write
        w('<?xml version="1.0" encoding="UTF-8"?>\n')
        w(f'<svg xmlns="http://www.w3.org/2000/svg" '
          f'xmlns:xlink="http://www.w3.org/1999/xlink" '
          f'width="{TOTAL_W}" height="{TOTAL_H}" viewBox="0 0 {TOTAL_W} {TOTAL_H}">\n')
        w(f'<rect width="{TOTAL_W}" height="{TOTAL_H}" fill="white"/>\n')

        for key, (x, y) in PANEL_POS.items():
            w(f'<g id="panel_{key}">')
            w(_svg_image(PANEL_FILES[key], x, y, PANEL_PX, PANEL_PX, embed, svg_dir))
            w("</g>\n")

        w('<g id="panel_E_bar">')
        w(_inline_plotly_svg(bar_svg, LEFT_W, 0, BAR_W_PX, BAR_H_PX))
        w("</g>\n")

        w('<g id="panel_letters">')
        for label, (x, y) in {**PANEL_POS, "E": (LEFT_W, 0)}.items():
            w(_svg_text(x + ox, y + oy, label, PANEL_LABEL_FONT_SIZE, bold=True))
        w(_svg_text(STRIP_X0 + ox, MAIN_H + oy, "F", PANEL_LABEL_FONT_SIZE, bold=True))
        w("</g>\n")

        w('<g id="panel_F_strip">\n')
        for i, group in enumerate(GROUPS_FIG6):
            r, c = divmod(i, STRIP_COLUMNS)
            x0 = STRIP_X0 + c * CELL_PX
            y0 = STRIP_Y0 + r * ROW_H
            sx, sy = x0 + HEADER_X_OFFSET, y0 + 14
            gid = re.sub(r"[^A-Za-z0-9_]", "_", group)

            w(f'<g id="F_{gid}">')
            w(f'<rect x="{sx}" y="{sy}" width="{HEADER_SWATCH_SIZE}" height="{HEADER_SWATCH_SIZE}" '
              f'fill="{group_colors[group]}" stroke="black" stroke-width="1"/>')
            w(_svg_text(sx + HEADER_SWATCH_SIZE + 14, sy - 4, group,
                        HEADER_NAME_FONT_SIZE, bold=True))
            st = stats_lookup.get(group)
            w(_svg_text(sx, sy + HEADER_SWATCH_SIZE + 14,
                        f"{st['Percent_inside']:.1f}% inside" if st is not None else "no data",
                        HEADER_STATS_FONT_SIZE))
            p = strip_panel_path(group)
            if os.path.exists(p):
                w(_svg_image(p, x0, y0 + STRIP_HEADER_HEIGHT, CELL_PX, CELL_PX, embed, svg_dir))
            w("</g>\n")
        w("</g>\n</svg>\n")

    print(f"Wrote pure-SVG figure: {out_path} ({TOTAL_W} x {TOTAL_H} units, "
          f"{os.path.getsize(out_path) / 1e6:.1f} MB)")
    return out_path


print("\n=========================================================")
print("ASSEMBLING PURE SVG")
print("=========================================================")

SVG_PATH = assemble_svg(os.path.join(OUT_DIR, FULL_NAME + ".svg"))


# ============================================================================
# RASTER COMPOSITE (PNG + TIFF) from the same fresh panels
# ============================================================================

def load_font(path, size):
    try:
        return ImageFont.truetype(path, size)
    except Exception:
        return ImageFont.load_default()


def _load_resized(path, size):
    im = Image.open(path).convert("RGB")
    if im.size != size:
        im = im.resize(size, Image.LANCZOS)
    return im


if EXPORT_RASTER_COMPOSITE:
    print("\n=========================================================")
    print("RASTER COMPOSITE (PNG + TIFF)")
    print("=========================================================")

    canvas = Image.new("RGB", (TOTAL_W, TOTAL_H), "white")
    draw = ImageDraw.Draw(canvas)
    big_font  = load_font(PANEL_LABEL_FONT_PATH, PANEL_LABEL_FONT_SIZE)
    name_font = load_font(HEADER_FONT_PATH, HEADER_NAME_FONT_SIZE)
    stat_font = load_font(HEADER_STATS_FONT_PATH, HEADER_STATS_FONT_SIZE)
    ox, oy = PANEL_LABEL_OFFSET

    for key, (x, y) in PANEL_POS.items():
        canvas.paste(_load_resized(PANEL_FILES[key], (PANEL_PX, PANEL_PX)), (x, y))
    canvas.paste(_load_resized(BAR_PNG_FILE, (BAR_W_PX, BAR_H_PX)), (LEFT_W, 0))

    for label, (x, y) in {**PANEL_POS, "E": (LEFT_W, 0)}.items():
        draw.text((x + ox, y + oy), label, fill="black", font=big_font)
    draw.text((STRIP_X0 + ox, MAIN_H + oy), "F", fill="black", font=big_font)

    for i, group in enumerate(GROUPS_FIG6):
        r, c = divmod(i, STRIP_COLUMNS)
        x0 = STRIP_X0 + c * CELL_PX
        y0 = STRIP_Y0 + r * ROW_H
        sx, sy = x0 + HEADER_X_OFFSET, y0 + 14
        draw.rectangle([sx, sy, sx + HEADER_SWATCH_SIZE, sy + HEADER_SWATCH_SIZE],
                       fill=group_colors[group], outline="black")
        draw.text((sx + HEADER_SWATCH_SIZE + 14, sy - 4), group, fill="black", font=name_font)
        st = stats_lookup.get(group)
        draw.text((sx, sy + HEADER_SWATCH_SIZE + 14),
                  f"{st['Percent_inside']:.1f}% inside" if st is not None else "no data",
                  fill="black", font=stat_font)
        p = strip_panel_path(group)
        if os.path.exists(p):
            canvas.paste(_load_resized(p, (CELL_PX, CELL_PX)), (x0, y0 + STRIP_HEADER_HEIGHT))

    png_out  = os.path.join(OUT_DIR, FULL_NAME + ".png")
    tiff_out = os.path.join(OUT_DIR, FULL_NAME + ".tiff")
    canvas.save(png_out, format="PNG")
    canvas.save(tiff_out, format="TIFF", dpi=(TIFF_DPI, TIFF_DPI), compression="tiff_lzw")
    print(f"Wrote {png_out} and {tiff_out} ({canvas.width} x {canvas.height})")


# ============================================================================
# OPTIONAL: Inkscape -> outlined-text SVG + PDF (fonts can't go missing)
# ============================================================================

if INKSCAPE_OUTLINED:
    ink = shutil.which("inkscape")
    if ink is None:
        print("\nInkscape not found on PATH — skipped outlined SVG/PDF "
              "(install Inkscape, or run it by hand on the SVG).")
    else:
        for ext, extra in (("pdf", []), ("svg", ["--export-plain-svg"])):
            out = os.path.join(OUT_DIR, f"{FULL_NAME}_outlined.{ext}" if ext == "svg"
            else f"{FULL_NAME}.pdf")
            cmd = [ink, SVG_PATH, f"--export-type={ext}", f"--export-filename={out}",
                   "--export-text-to-path"] + extra
            print("Running:", " ".join(cmd))
            try:
                subprocess.run(cmd, check=True)
            except Exception as e:
                print(f"  Inkscape export to {ext} failed: {e}")


print("\n=========================================================")
print("DONE")
print("=========================================================")
print(f"Pure SVG        : {os.path.abspath(SVG_PATH)}")
print(f"Panels (fresh)  : {os.path.abspath(PANEL_DIR)}/")
print(f"Other exports   : {os.path.abspath(OUT_DIR)}/")
print(f"Stats CSV       : {STATS_CSV}")