#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure4_final.py   (v2: PATTERN version)
Final Figure 4 (7 rows x 7 columns) + population summary + supplementary figures
=================================================================================

STRATEGY (nothing is shuffled, the - blocks are never reversed / reverse-complemented)
  A "block" = one unbroken run of one annotation class on one strand; its A,T,G,C counts are kept in
  reference coordinates.  Blocks are grouped in 10-nt LENGTH BINS (CDS: about 40 bins for 20-420 nt).
  Everything below is a PATTERN ALONG THE LENGTH BINS (which bins are high, which are low):

  1. LENGTH PATTERN   % of blocks in every length bin, + vs -.
                      D = KS distance between the two length distributions, r = Spearman correlation of
                      the two bin-frequency patterns.
  2. NUCLEOTIDE PATTERN  portion of a base in every length bin (bases pooled over all blocks of the bin):
                      A in +  vs  T in -,   T in + vs A in -,   G in + vs C in -,   C in + vs G in -.
                      r = Spearman correlation of the two patterns across the bins (do they go up and down
                      together?),  Delta = base-weighted mean |difference| of the two patterns (pp).
  3. CROSS TEST       the eight patterns A+, T+, C+, G+, A-, T-, C-, G-: all 8 x 8 correlations (diagonal = 1 =
                      maximum, "A+ with A+") and all 8 x 8 distances (diagonal = 0).
                      Baseline "self": one pattern built from odd-numbered blocks vs the same pattern built from
                      even-numbered blocks (no random shuffling) = what an identical pattern gives with this sampling.
  4. NEAREST          is the RC partner (A+ -> T-, T+ -> A-, ...) the nearest of the other patterns?
                      rank 1 = nearest (highest r / smallest Delta), among the 7 others.
  5. RATIOS           per length bin from reference counts: A+/T+ (+ only), A-/T- (- only) and the pooled ratio
                      (A+ + A-)/(T+ + T-); same for G/C.  Deviation = base-weighted mean |R-1|.
                      Silencing = mean(deviation of + only, - only) / deviation of the pooled ratio.

  "x chance" = D / (D expected for two samples of the same sizes from ONE distribution, 0.87*sqrt((n+m)/(n*m));
  continuous-data approximation, ties make it slightly smaller).

MODES
  --mode main        Human_T2T chromosome 1 (first valid chromosome, or --chromosome ID)
                       -> Figure_4_<organism>_<chr>.(pdf|svg|tiff|png)   7 x 7 figure, statistics in every panel
                       -> Supp_Fig_chr1_pattern_matrix.(pdf|png)         8 x 8 r and Delta matrices per class
                       -> Figure4_stats_<chr>.csv, _pairs_, _ratios_, _cross_matrices_, Figure4_report_<chr>.txt
  --mode population  every chromosome of every organism folder (streamed, parallel, resumable)
                       -> per_chromosome_similarity.csv, summary_by_class.csv, Population_report.txt
                       -> Supp_Fig_S1_population_7col   (7 columns like Figure 4: length, 4 RC pairs, A/T, G/C;
                                                          3 statistics per column; one box = all chromosomes of a class)
                       -> Supp_Fig_S2_population_cross_matrix  (median 8 x 8 pattern matrices per class)
                       -> with --extended also S3-S6 (length, pairs, cross/nearest, ratios in detail)
  --mode summarise   rebuild the population tables / report / figures from finished organisms
  --mode all         main, then population

LABELS (Figure 4 panel titles)   --labels short (default) | full
  short : minimal statistics for the main figure  -  length: D (x chance), N+/N-;  RC pairs: Delta (pp) and
          "RC nearest" (tick = the RC partner is the nearest of the other 7 patterns by Delta, otherwise its rank);
          ratios: strand-only deviations + | - and the pooled deviation after the arrow.
  full  : every statistic (r, ranks by r and Delta, silencing factor) - use it for checking / supplement.
  All statistics are always written to the CSV files and to the Excel workbooks (Supp_Tables_*.xlsx).

FIGURE 4 COLUMNS (all panels of one row share the length axis)
  1  length distribution (% of blocks per 10-nt bin, - left | + right)      D, x chance, r, N+/N-
  2-5 nucleotide pattern per length bin  T(-)|A(+), A(-)|T(+), C(-)|G(+), G(-)|C(+)   r, Delta, rank of RC partner
  6  A/T ratio per length bin: + only, - only, pooled (grey band 95% CI)      deviations and silencing
  7  G/C ratio per length bin: same

Usage
  python Figure4_final.py --mode main
  python Figure4_final.py --mode population --workers 4
  python Figure4_final.py --mode summarise
  python Figure4_final.py --mode main --chromosome NC_060925.1 --organism Human_T2T
"""

from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import argparse
import gc
import glob
import gzip
import hashlib
import json
import multiprocessing as mp
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

VERSION = "2.2"

# =============================================================================
# 1. INPUT / OUTPUT
# =============================================================================

INPUT_ROOT = str(config.MASKED_FASTAS)   # Fasta_GFF_inputs/Masked_Fastas/<organism>/
OUTPUT_DIR = str(config.output_dir(4))   # Figure_4/outputs
ORGANISM = "Human_T2T"
MAIN_CHROMOSOME: str | None = None          # None = first valid chromosome of ORGANISM (chromosome 1)
MIN_CHROMOSOME_BP = 5_000_000
MAX_FEATURE_LENGTH = 2500
CHUNK_MB = 16
LINE_LIMIT = 1 << 24
WORKERS = min(4, os.cpu_count() or 1)

# =============================================================================
# 2. WINDOWS / BINS (identical to Figure 4 v2.3)
# =============================================================================

# ---- LENGTH RANGE OF EVERY CLASS (nt): set it here, by hand -------------------------------
# One entry per class. These are the analysed windows and (unless DISPLAY_Y_RANGE says otherwise)
# the length axis of the figure.
ROW_Y_MIN = {"LINE": 20, "SINE": 20, "LTR": 20, "DNA": 20, "CDS": 20, "INTRON": 20, "UTR": 20}
ROW_Y_MAX = {"LINE": 2000, "SINE": 400, "LTR": 1600, "DNA": 1000, "CDS": 820, "INTRON": 2500, "UTR": 2000}
# ONE range for ALL classes (overrides the two dictionaries above), e.g. (20, 1000). None = use the dictionaries.
SAME_RANGE_FOR_ALL = None
# Figure only: show another length axis for a class without changing the analysis, e.g. {"LINE": (20, 1000)}.
DISPLAY_Y_RANGE: dict = {}


def class_window(feature):
    """(lowest, highest) block length in nt analysed for this class."""
    if SAME_RANGE_FOR_ALL is not None:
        return int(SAME_RANGE_FOR_ALL[0]), int(SAME_RANGE_FOR_ALL[1])
    return int(ROW_Y_MIN.get(feature, 20)), int(ROW_Y_MAX.get(feature, MAX_FEATURE_LENGTH or 2500))


BIN_WIDTH_NT = 10
MIN_ELEMENTS_PER_STRAND_BIN = 10
TRIM_LENGTH_RANGE = True             # True: cut the sparse tail automatically; False: use exactly the range set above (--no-trim)
MAX_GAP_BINS = 3
MIN_BLOCKS_PER_STRAND = 50           # fewer blocks -> class reported as insufficient
MIN_USABLE_BINS = 8                  # fewer usable length bins -> class reported as insufficient (patterns need bins)
KS_MEAN = 0.8687                     # mean of the Kolmogorov distribution = sqrt(pi/2) * ln 2
SQRT_2_OVER_PI = float(np.sqrt(2.0 / np.pi))
CI_Z = 1.96

# =============================================================================
# 3. FIGURE GEOMETRY (Nucleic Acids Research: 17.5 cm = double column)
# =============================================================================

FIG_WIDTH_CM, FIG_HEIGHT_CM = 17.5, 23.0
MARGIN_LEFT, MARGIN_RIGHT, MARGIN_TOP, MARGIN_BOTTOM = 0.115, 0.985, 0.992, 0.060
WSPACE, HSPACE, HEADER_RATIO = 0.22, 0.85, 0.80
COL_WIDTH_RATIOS = [1.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
HEADER_GROUP_LABEL_CM, HEADER_GROUP_RULE_CM = 0.62, 0.95
HEADER_COL_TITLE_CM, HEADER_COL_SUBTITLE_CM = 1.10, 1.45
TOP_LEGEND_DROP = -0.004
FS_TICK, FS_LABEL, FS_STAT, FS_LEGEND = 6.0, 6.5, 5.5, 6.0
FS_COL_TITLE, FS_SUB, FS_GROUP, FS_ROW = 7.5, 5.8, 7.0, 8.5
FONT_STACK = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]
CI_ALPHA = 0.16

# =============================================================================
# 4. FEATURE ENCODING (unchanged)
# =============================================================================

FEATURE_TYPES = {
    "LINE":   {"+": "*YUI", "-": "+yui"},
    "SINE":   {"+": "OPJK", "-": "opjk"},
    "LTR":    {"+": "LM@#", "-": "lm$%"},
    "DNA":    {"+": "1234", "-": "5678"},
    "CDS":    {"+": "ZXVB", "-": "zxvb"},
    "INTRON": {"+": "QWER", "-": "qwer"},
    "UTR":    {"+": "SDFH", "-": "sdfh"},
}
FEATURES = list(FEATURE_TYPES)
COLORS = {"LINE": ("#1f77b4", "#9ecae1"), "SINE": ("#d62728", "#f7b6b2"), "LTR": ("#2ca02c", "#98df8a"),
          "DNA": ("#9467bd", "#c5b0d5"), "CDS": ("#ff7f0e", "#ffbb78"), "INTRON": ("#17becf", "#9edae5"),
          "UTR": ("#8c564b", "#c49c94")}

# counts columns are A, T, G, C.  Main pairs: (index in +, index in -, key, figure column title)
PAIRS = [(0, 1, "AplusTminus", "T(\u2212) | A(+)", "A(+) vs T(\u2212)"),
         (1, 0, "TplusAminus", "A(\u2212) | T(+)", "T(+) vs A(\u2212)"),
         (2, 3, "GplusCminus", "C(\u2212) | G(+)", "G(+) vs C(\u2212)"),
         (3, 2, "CplusGminus", "G(\u2212) | C(+)", "C(+) vs G(\u2212)")]
PAIR_KEYS = [p[2] for p in PAIRS]
PAIR_COLORS = {"AplusTminus": "#0072B2", "TplusAminus": "#56B4E9",
               "GplusCminus": "#D55E00", "CplusGminus": "#E69F00"}
# cross test: the eight histograms in this order
DISTS = [("A", 0, "+"), ("T", 1, "+"), ("C", 3, "+"), ("G", 2, "+"),
         ("A", 0, "-"), ("T", 1, "-"), ("C", 3, "-"), ("G", 2, "-")]
DIST_LABELS = [f"{b}{s}" for b, _, s in DISTS]
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def pretty(lab):
    return lab.replace("-", "\u2212")


FEAT_CODE = np.full(256, -1, dtype=np.int8)
BASE_CODE = np.full(256, -1, dtype=np.int8)
for _fi, (_n, _info) in enumerate(FEATURE_TYPES.items()):
    for _si, _s in enumerate(("+", "-")):
        _chars = _info[_s]
        if len(_chars) != 4 or len(set(_chars)) != 4:
            raise ValueError(f"{_n} {_s}: need 4 distinct characters (A,T,G,C)")
        for _bi, _ch in enumerate(_chars):
            if FEAT_CODE[ord(_ch)] != -1:
                raise ValueError(f"Encoding collision on {_ch!r}")
            FEAT_CODE[ord(_ch)] = _fi * 2 + _si
            BASE_CODE[ord(_ch)] = _bi


# =============================================================================
# 5. FASTA DISCOVERY + STREAMING BLOCK EXTRACTION
# =============================================================================

def normalize_text(t):
    return re.sub(r"\s+", " ", str(t).lower()).strip()


def safe_name(t):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(t))


def find_fasta_files(genome_dir):
    out = []
    for root, _, files in os.walk(genome_dir):
        for f in files:
            if f.lower().endswith((".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")):
                out.append(os.path.join(root, f))
    return sorted(out)


def get_organism_dirs(root, targets):
    d = [p for p in Path(root).iterdir() if p.is_dir()]
    if targets:
        d = [p for p in d if p.name in set(targets)]
    return sorted(str(p) for p in d)


def is_valid_chromosome_name(seq_id, description, seq_len):
    text = normalize_text(f"{seq_id} {description}")
    if seq_len < MIN_CHROMOSOME_BP:
        return False
    if any(w in text for w in ["scaffold", "contig", "unlocalized", "unplaced", "patch",
                               "plasmid", "mitochond", "chloroplast"]):
        return False
    if "chromosome un" in text or re.search(r"\bchr\s*un\b", text):
        return False
    if re.search(r"\bchr[\s._-]*[0-9xyzwm]+\b", text):
        return True
    if re.match(r"^(NC_|CM\d|CP\d)", seq_id.upper()):
        return True
    if "chromosome" in text:
        return True
    if "complete sequence" in text or "complete genome" in text:
        return True
    return False


class BlockCollector:
    """Builds the annotation blocks of one chromosome chunk by chunk (reference coordinates)."""

    def __init__(self):
        self.parts = {(fi, si): [] for fi in range(len(FEATURES)) for si in (0, 1)}
        self.open_f = -1
        self.open_c = np.zeros(4, np.int64)

    def _store(self, run_f, counts):
        keep = run_f >= 0
        if not keep.any():
            return
        rf, c = run_f[keep], counts[keep]
        if MAX_FEATURE_LENGTH is not None:
            ok = c.sum(1) <= MAX_FEATURE_LENGTH
            rf, c = rf[ok], c[ok]
        for code in np.unique(rf):
            self.parts[(int(code) // 2, int(code) % 2)].append(c[rf == code].astype(np.int32))

    def feed(self, arr):
        if arr.size == 0:
            return
        f = FEAT_CODE[arr]
        starts = np.concatenate(([0], np.flatnonzero(f[1:] != f[:-1]) + 1))
        run_f = f[starts]
        base = BASE_CODE[arr]
        counts = np.empty((len(starts), 4), np.int64)
        for k in range(4):
            counts[:, k] = np.add.reduceat(base == k, starts, dtype=np.int64)
        del f, base
        if self.open_f == run_f[0]:
            counts[0] += self.open_c                       # run continues across the chunk border
        else:
            self._store(np.array([self.open_f]), self.open_c[None, :])
        self._store(run_f[:-1], counts[:-1])
        self.open_f, self.open_c = int(run_f[-1]), counts[-1].copy()

    def result(self):
        self._store(np.array([self.open_f]), self.open_c[None, :])
        out = {}
        for fi, name in enumerate(FEATURES):
            pair = []
            for si in (0, 1):
                p = self.parts[(fi, si)]
                pair.append(np.concatenate(p) if p else np.zeros((0, 4), np.int32))
            out[name] = tuple(pair)
        self.parts = None
        return out


def iter_chromosomes(path, wanted, cov, seen=None):
    """Yield (seq_id, header, blocks, length) for valid chromosomes of one FASTA file.
    If `wanted` is a set of ids, other chromosomes are skipped without being parsed."""
    seen = set() if seen is None else seen
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    chunk = max(1, int(CHUNK_MB * 1024 * 1024))
    col, buf, total, cur_id, cur_header = None, bytearray(), 0, None, None
    cont = False

    def close(col, buf, total, cur_id, cur_header):
        if col is None:
            return None
        if buf:
            col.feed(np.frombuffer(bytes(buf), np.uint8))
            buf.clear()
        return (cur_id, cur_header, col.result(), total)

    def accept(rec):
        if rec is None:
            return False
        if rec[3] >= MIN_CHROMOSOME_BP:
            return True
        cov["records_too_short"] += 1
        return False

    with opener(path, "rb") as fh:
        while True:
            raw = fh.readline(LINE_LIMIT)
            if not raw:
                break
            was_cont, cont = cont, not raw.endswith(b"\n")
            if not was_cont and raw.startswith(b">"):
                rec = close(col, buf, total, cur_id, cur_header)
                col = None
                if accept(rec):
                    yield rec
                rec = None
                header = raw[1:].decode("ascii", "replace").strip()
                seq_id = header.split()[0] if header.split() else header
                cov["records_total"] += 1
                if not is_valid_chromosome_name(seq_id, header, MIN_CHROMOSOME_BP):
                    cov["records_rejected_name"] += 1
                elif seq_id in seen:
                    cov["duplicate_ids_skipped"] += 1
                elif wanted is not None and seq_id not in wanted:
                    cov["records_not_wanted"] += 1
                else:
                    seen.add(seq_id)
                    col, buf, total, cur_id, cur_header = BlockCollector(), bytearray(), 0, seq_id, header
                continue
            if col is not None:
                line = raw.strip()
                buf += line
                total += len(line)
                if len(buf) >= chunk:
                    col.feed(np.frombuffer(bytes(buf), np.uint8))
                    buf.clear()
    rec = close(col, buf, total, cur_id, cur_header)
    if accept(rec):
        yield rec


# =============================================================================
# 6. BLOCK SELECTION (same window / bins / tail trimming as Figure 4)
# =============================================================================

def bin_geometry(feature):
    lo, hi = class_window(feature)
    w = max(1, int(BIN_WIDTH_NT))
    nb = (hi - lo + 1) // w
    return dict(lo=lo, hi_eff=lo + nb * w - 1, w=w, nb=nb)


def select_blocks(feature, cp, cm):
    g = bin_geometry(feature)
    lo, hi_eff, w, nb = g["lo"], g["hi_eff"], g["w"], g["nb"]
    Lp, Lm = cp.sum(1), cm.sum(1)
    sp_, sm_ = (Lp >= lo) & (Lp <= hi_eff), (Lm >= lo) & (Lm <= hi_eff)
    cp, cm, Lp, Lm = cp[sp_], cm[sm_], Lp[sp_], Lm[sm_]
    ip, im = ((Lp - lo) // w).astype(int), ((Lm - lo) // w).astype(int)
    cut = nb - 1
    ok0 = ((np.bincount(ip, minlength=nb) >= MIN_ELEMENTS_PER_STRAND_BIN)
           & (np.bincount(im, minlength=nb) >= MIN_ELEMENTS_PER_STRAND_BIN))
    if TRIM_LENGTH_RANGE and ok0.any():
        usable = np.flatnonzero(ok0)
        cut = int(usable[0])
        for b in usable[1:]:
            if b - cut - 1 <= MAX_GAP_BINS:
                cut = int(b)
            else:
                break
    kp, km = ip <= cut, im <= cut
    return cp[kp], cm[km], ip[kp], im[km], nb, cut, g


# =============================================================================
# 7. STATISTICS (patterns along the length bins, cross test, nearest pattern, ratios)
# =============================================================================

def ks_sorted(xs, ys):
    """Two-sample KS distance D for already sorted samples (exact, ties handled)."""
    if len(xs) == 0 or len(ys) == 0:
        return np.nan
    allv = np.concatenate([xs, ys])
    f1 = np.searchsorted(xs, allv, side="right") / len(xs)
    f2 = np.searchsorted(ys, allv, side="right") / len(ys)
    return float(np.max(np.abs(f1 - f2)))


def ks_d(x, y):
    return ks_sorted(np.sort(x), np.sort(y))


def ks_expected(n, m):
    return KS_MEAN * np.sqrt((n + m) / (n * m)) if n > 0 and m > 0 else np.nan


def _bc(idx, weights, nb):
    return np.bincount(idx, weights=weights, minlength=nb)


def ratio_by_bin(a, b, idx, nb):
    """Pooled ratio sum(a)/sum(b) per bin with cluster-robust SE (delta method, block = unit)."""
    n = np.bincount(idx, minlength=nb).astype(float)
    a, b = a.astype(float), b.astype(float)
    A, B = _bc(idx, a, nb), _bc(idx, b, nb)
    with np.errstate(divide="ignore", invalid="ignore"):
        R = A / B
        e = a - R[idx] * b
        V = _bc(idx, e ** 2, nb)
        se = np.sqrt(V / B ** 2 * (n / (n - 1.0)))
    return R, se


def _wmean(values, weights, mask):
    v = np.where(mask, values, np.nan)
    m = np.isfinite(v) & np.isfinite(weights) & (weights > 0)
    return float(np.average(v[m], weights=weights[m])) if m.any() else np.nan


def rank_corr(x, y):
    """Spearman correlation of two patterns (NaN if constant or too short)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan
    rx, ry = sps.rankdata(x[m]), sps.rankdata(y[m])
    if rx.std() == 0 or ry.std() == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def rank_corr_matrix(prof):
    ranks = np.vstack([sps.rankdata(v) for v in prof])
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.corrcoef(ranks)


def dist_matrix(prof, w):
    """Base-weighted mean absolute difference (pp) between all pairs of patterns."""
    D = np.zeros((len(prof), len(prof)))
    for i in range(len(prof)):
        D[i] = (np.abs(prof[i][None, :] - prof) * w[None, :]).sum(1) / w.sum()
    return D


def half_profiles(C, idx, nb):
    """Pattern built from odd-numbered blocks and from even-numbered blocks (genomic order, no shuffling)."""
    par = np.arange(len(C)) % 2
    out = []
    for h in (0, 1):
        m = par == h
        S = np.column_stack([_bc(idx[m], C[m, k], nb) for k in range(4)])
        with np.errstate(invalid="ignore", divide="ignore"):
            out.append(S / S.sum(1)[:, None] * 100)
    return out


def partner_index(i):
    base, _, strand = DISTS[i]
    return DIST_LABELS.index(f"{COMP[base]}{'-' if strand == '+' else '+'}")


def _rank(row, i, p, cand, higher_is_better):
    x = row[p]
    if not np.isfinite(x):
        return np.nan
    v = row[cand]
    return 1 + int(((v > x) if higher_is_better else (v < x)).sum())


BASES_ATGC = "ATGC"
DKEYS = ["Ap", "Tp", "Cp", "Gp", "Am", "Tm", "Cm", "Gm"]      # column-name keys of the eight patterns


def analyse_feature(feature, cp, cm, keep=False):
    """All statistics of one class on one chromosome. keep=True also returns the arrays for plotting."""
    cp, cm, ip, im, nb, cut, g = select_blocks(feature, cp, cm)
    n_p, n_m = len(cp), len(cm)
    res = {"feature": feature, "status": "insufficient", "n_p": n_p, "n_m": n_m,
           "pairs": [], "ratios": [], "cross": [], "arr": None,
           "wide": {"Feature": feature, "Status": "insufficient", "N_plus": n_p, "N_minus": n_m}}
    if n_p < MIN_BLOCKS_PER_STRAND or n_m < MIN_BLOCKS_PER_STRAND:
        return res
    cnt_p, cnt_m = np.bincount(ip, minlength=nb), np.bincount(im, minlength=nb)
    used = np.flatnonzero((cnt_p >= MIN_ELEMENTS_PER_STRAND_BIN) & (cnt_m >= MIN_ELEMENTS_PER_STRAND_BIN))
    if len(used) < MIN_USABLE_BINS:
        return res
    res["status"] = "ok"
    wide = {"Feature": feature, "Status": "ok", "N_plus": n_p, "N_minus": n_m, "Bins_used": int(len(used))}

    # ---- 1. LENGTH pattern: % of blocks per length bin -----------------------------------
    Lp, Lm = cp.sum(1).astype(float), cm.sum(1).astype(float)
    d_len, e_len = ks_d(Lp, Lm), ks_expected(n_p, n_m)
    fp, fm = cnt_p[:cut + 1] / n_p * 100, cnt_m[:cut + 1] / n_m * 100
    wide.update({"Len_balance_pct": (1 - abs(n_p - n_m) / (n_p + n_m)) * 100, "Len_KS_D": d_len,
                 "Len_KS_D_expected": e_len, "Len_D_over_expected": d_len / e_len,
                 "Len_r": rank_corr(fp, fm)})

    # ---- 2. NUCLEOTIDE patterns: pooled portion of each base in every length bin ----------------
    Sp = np.column_stack([_bc(ip, cp[:, k], nb) for k in range(4)])
    Sm = np.column_stack([_bc(im, cm[:, k], nb) for k in range(4)])
    Np, Nm = Sp.sum(1), Sm.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        Pp, Pm = Sp / Np[:, None] * 100, Sm / Nm[:, None] * 100          # percent of bases in the bin
    prof = np.array([(Pp if s == "+" else Pm)[used, k] for _, k, s in DISTS])     # 8 patterns
    w = (Np + Nm)[used]

    # ---- 3. CROSS TEST: 8 x 8 correlations and distances, self baseline -------------------------
    R = rank_corr_matrix(prof)
    D = dist_matrix(prof, w)
    for i_ in range(8):
        for j_ in range(i_ + 1, 8):
            wide[f"R_{DKEYS[i_]}_{DKEYS[j_]}"] = R[i_, j_]
            wide[f"D_{DKEYS[i_]}_{DKEYS[j_]}"] = D[i_, j_]
    Hp, Hm = half_profiles(cp, ip, nb), half_profiles(cm, im, nb)
    r_self, d_self = np.full(8, np.nan), np.full(8, np.nan)
    for i, (_, k, s) in enumerate(DISTS):
        h0, h1 = (Hp if s == "+" else Hm)[0][used, k], (Hp if s == "+" else Hm)[1][used, k]
        r_self[i] = rank_corr(h0, h1)
        mk = np.isfinite(h0) & np.isfinite(h1)
        d_self[i] = _wmean(np.abs(h0 - h1), w, mk)

    # ---- 4. NEAREST pattern: rank of the RC partner ------------------------------------------------
    rk = {"r7": np.full(8, np.nan), "d7": np.full(8, np.nan), "r4": np.full(8, np.nan), "d4": np.full(8, np.nan)}
    best_r, best_d, part_r, part_d = np.full(8, np.nan), np.full(8, np.nan), np.full(8, np.nan), np.full(8, np.nan)
    for i in range(8):
        p = partner_index(i)
        others = [j for j in range(8) if j != i]
        other_strand = [4, 5, 6, 7] if DISTS[i][2] == "+" else [0, 1, 2, 3]
        non = [j for j in range(8) if j not in (i, p)]
        rk["r7"][i] = _rank(R[i], i, p, others, True)
        rk["d7"][i] = _rank(D[i], i, p, others, False)
        rk["r4"][i] = _rank(R[i], i, p, other_strand, True)
        rk["d4"][i] = _rank(D[i], i, p, other_strand, False)
        part_r[i], part_d[i] = R[i, p], D[i, p]
        best_r[i] = np.nanmax(R[i, non]) if np.isfinite(R[i, non]).any() else np.nan
        best_d[i] = np.nanmin(D[i, non])
        res["cross"].append({
            "Feature": feature, "Row": DIST_LABELS[i], "RC_partner": DIST_LABELS[p],
            "r_partner": R[i, p], "r_self_split": r_self[i], "r_best_nonpartner": best_r[i],
            "Delta_partner": D[i, p], "Delta_self_split": d_self[i], "Delta_best_nonpartner": best_d[i],
            "Nearest_by_r": DIST_LABELS[others[int(np.nanargmax(R[i, others]))]]
            if np.isfinite(R[i, others]).any() else "",
            "Nearest_by_Delta": DIST_LABELS[others[int(np.argmin(D[i, others]))]],
            "Rank_r_among_7": rk["r7"][i], "Rank_Delta_among_7": rk["d7"][i],
            "Rank_r_among_4_other_strand": rk["r4"][i], "Rank_Delta_among_4_other_strand": rk["d4"][i]})
        wide[f"RankR7_{DIST_LABELS[i]}"] = rk["r7"][i]
        wide[f"RankD7_{DIST_LABELS[i]}"] = rk["d7"][i]
    wide.update({"Cross_r_nearest_of7": int(np.nansum(rk["r7"] == 1)), "Cross_d_nearest_of7": int(np.nansum(rk["d7"] == 1)),
                 "Cross_r_nearest_of4": int(np.nansum(rk["r4"] == 1)), "Cross_d_nearest_of4": int(np.nansum(rk["d4"] == 1)),
                 "Cross_r_partner_median": float(np.nanmedian(part_r)), "Cross_r_self_median": float(np.nanmedian(r_self)),
                 "Cross_r_nonpartner_median": float(np.nanmedian(best_r)),
                 "Cross_d_partner_median": float(np.nanmedian(part_d)), "Cross_d_self_median": float(np.nanmedian(d_self)),
                 "Cross_d_nonpartner_median": float(np.nanmedian(best_d))})

    # ---- the four target pairs (A+ vs T-, T+ vs A-, G+ vs C-, C+ vs G-) -----------------------------
    for a, c, key, _, name in PAIRS:
        i_p, i_m = DIST_LABELS.index(f"{BASES_ATGC[a]}+"), DIST_LABELS.index(f"{BASES_ATGC[c]}-")
        rr = max(rk["r7"][i_p], rk["r7"][i_m]) if np.isfinite([rk["r7"][i_p], rk["r7"][i_m]]).all() else np.nan
        rd = max(rk["d7"][i_p], rk["d7"][i_m]) if np.isfinite([rk["d7"][i_p], rk["d7"][i_m]]).all() else np.nan
        row = {"Feature": feature, "Pair": name, "Bins_used": int(len(used)),
               "Mean_plus_pct": 100 * Sp[used, a].sum() / Np[used].sum(),
               "Mean_minus_pct": 100 * Sm[used, c].sum() / Nm[used].sum(),
               "r": R[i_p, i_m], "Delta_pp": D[i_p, i_m],
               "r_self": float(np.nanmean([r_self[i_p], r_self[i_m]])),
               "Delta_self_pp": float(np.nanmean([d_self[i_p], d_self[i_m]])),
               "Rank_r_among_7": rr, "Rank_Delta_among_7": rd}
        res["pairs"].append(row)
        wide.update({f"{key}_r": row["r"], f"{key}_delta": row["Delta_pp"], f"{key}_r_self": row["r_self"],
                     f"{key}_delta_self": row["Delta_self_pp"], f"{key}_rank_r": rr, f"{key}_rank_d": rd,
                     f"{key}_mean_plus": row["Mean_plus_pct"], f"{key}_mean_minus": row["Mean_minus_pct"]})

    # ---- 5. RATIOS per length bin: + only, - only, pooled -----------------------------------------------
    idx_all, ref_all = np.concatenate([ip, im]), np.vstack([cp, cm])
    mask_used = np.isin(np.arange(nb), used)
    ratio_arr = {}
    for key, (a, b) in (("AT", (0, 1)), ("GC", (2, 3))):
        with np.errstate(invalid="ignore", divide="ignore"):
            Rp, Rm = Sp[:, a] / Sp[:, b], Sm[:, a] / Sm[:, b]
        Rt, se_t = ratio_by_bin(ref_all[:, a], ref_all[:, b], idx_all, nb)
        wts = (Sp[:, a] + Sp[:, b] + Sm[:, a] + Sm[:, b])
        dev_p, dev_m = _wmean(np.abs(Rp - 1), wts, mask_used), _wmean(np.abs(Rm - 1), wts, mask_used)
        dev_t = _wmean(np.abs(Rt - 1), wts, mask_used)
        noise = _wmean(SQRT_2_OVER_PI * se_t, wts, mask_used)
        sil = 0.5 * (dev_p + dev_m) / max(dev_t, 1e-9)
        r_opp = rank_corr((Rp - 1)[used], (Rm - 1)[used])
        row = {"Feature": feature, "Ratio": "A/T" if key == "AT" else "G/C", "Dev_plus": dev_p, "Dev_minus": dev_m,
               "Dev_pooled": dev_t, "Pooled_noise": noise, "Silencing": sil, "r_opposite": r_opp,
               "Rplus_median": float(np.nanmedian(Rp[used])), "Rminus_median": float(np.nanmedian(Rm[used]))}
        res["ratios"].append(row)
        wide.update({f"{key}_dev_plus": dev_p, f"{key}_dev_minus": dev_m, f"{key}_dev_pooled": dev_t,
                     f"{key}_noise": noise, f"{key}_silencing": sil, f"{key}_r_opp": r_opp,
                     f"{key}_Rplus_median": row["Rplus_median"], f"{key}_Rminus_median": row["Rminus_median"]})
        ratio_arr[key] = (Rp[used], Rm[used], Rt[used], se_t[used])
    res["wide"] = wide
    if keep:
        edges = g["lo"] + g["w"] * np.arange(cut + 2, dtype=float)
        res["arr"] = {"edges": edges, "len_p": fp, "len_m": fm, "used": used,
                      "ycent": g["lo"] + g["w"] * (used + 0.5), "bin_w": g["w"],
                      "Pp": Pp[used], "Pm": Pm[used], "R": R, "D": D, "r_self": r_self, "d_self": d_self,
                      "rank": rk, "ratio": ratio_arr}
    return res


# =============================================================================
# 8. FIGURE 4 (7 rows x 7 columns)
# =============================================================================

def apply_style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": FONT_STACK, "font.size": FS_TICK,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                         "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.pad": 1.5,
                         "ytick.major.pad": 1.5, "pdf.fonttype": 42, "ps.fonttype": 42,
                         "svg.fonttype": "path", "axes.unicode_minus": True})


LABELS_MODE = "short"        # "short" = minimal statistics in the panels; "full" = all statistics
WRITE_EXCEL = True           # write the supplementary tables as .xlsx (needs openpyxl)


def label_length(W, d):
    base = f"D {W['Len_KS_D']:.3f} (\u00d7{fmt_ratio(W['Len_D_over_expected'])})"
    if LABELS_MODE == "full":
        base += f" | r {W['Len_r']:.2f}"
    return base + f"\n+{d['n_p']:,} / \u2212{d['n_m']:,}"


def label_pair(W, key):
    rd = W[f"{key}_rank_d"]
    if LABELS_MODE == "full":
        return (f"r {W[f'{key}_r']:.2f} | \u0394 {W[f'{key}_delta']:.2f}\n"
                f"RC rank {fmt_rank(W[f'{key}_rank_r'])}|{fmt_rank(rd)} of 7")
    nearest = "RC nearest $\\checkmark$" if (rd is not None and np.isfinite(rd) and rd == 1) else f"RC rank {fmt_rank(rd)}"
    return f"\u0394 {W[f'{key}_delta']:.2f} pp\n{nearest}"


def label_ratio(W, key):
    dp, dm, dt = W[f"{key}_dev_plus"], W[f"{key}_dev_minus"], W[f"{key}_dev_pooled"]
    if LABELS_MODE == "full":
        return (f"dev +{dp:.2f} | \u2212{dm:.2f}\n"
                f"pooled {dt:.3f} (\u00f7{fmt_sil(W[f'{key}_silencing'])})")
    return f"{dp:.2f} | {dm:.2f} \u2192 {dt:.3f}"


def fmt_ratio(r):
    return f"{r:.1f}" if r < 10 else f"{r:.0f}"


def fmt_rank(v):
    return "\u2013" if v is None or not np.isfinite(v) else f"{int(v)}"


def fmt_sil(s):
    return ">99" if s > 99 else (f"{s:.0f}" if s >= 10 else f"{s:.1f}")


def draw_figure4(res, chrom, out_dir, base, footnote=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator

    apply_style()
    n_rows, n_cols = len(FEATURES), 7
    fig = plt.figure(figsize=(FIG_WIDTH_CM / 2.54, FIG_HEIGHT_CM / 2.54), facecolor="white")
    gs = fig.add_gridspec(n_rows + 1, n_cols, height_ratios=[HEADER_RATIO] + [1] * n_rows,
                          width_ratios=COL_WIDTH_RATIOS, left=MARGIN_LEFT, right=MARGIN_RIGHT,
                          top=MARGIN_TOP, bottom=MARGIN_BOTTOM, wspace=WSPACE, hspace=HSPACE)
    specs = [("Blocks per length bin", "% of blocks, \u2212 left | + right"),
             (PAIRS[0][3], "portion per bin"), (PAIRS[1][3], "portion per bin"),
             (PAIRS[2][3], "portion per bin"), (PAIRS[3][3], "portion per bin"),
             ("A/T ratio", "per length bin"), ("G/C ratio", "per length bin")]
    xlabels = ["% of blocks per bin\n\u2190 minus | plus \u2192"] + \
              ["% of bases in bin\n\u2190 minus | plus \u2192"] * 4 + ["A/T ratio", "G/C ratio"]

    hax = []
    for j in range(n_cols):
        ax = fig.add_subplot(gs[0, j])
        ax.axis("off")
        hax.append(ax)
    y_top = hax[0].get_position().y1
    y_of = lambda cm: y_top - cm / FIG_HEIGHT_CM
    for j, (title, sub) in enumerate(specs):
        pos = hax[j].get_position()
        xc = 0.5 * (pos.x0 + pos.x1)
        fig.text(xc, y_of(HEADER_COL_TITLE_CM), title, ha="center", va="top", fontsize=FS_COL_TITLE, fontweight="bold")
        fig.text(xc, y_of(HEADER_COL_SUBTITLE_CM), sub, ha="center", va="top", fontsize=FS_SUB)
    groups = [(0, 0, "Length pattern"),
              (1, 4, "Nucleotide pattern across length bins: RC pairs"),
              (5, 6, "Base ratios: deviation, silencing")]
    for c0, c1, label in groups:
        x0, x1 = hax[c0].get_position().x0, hax[c1].get_position().x1
        pad = 0.012 * (x1 - x0) if c1 > c0 else 0.04 * (x1 - x0)
        fig.text(0.5 * (x0 + x1), y_of(HEADER_GROUP_LABEL_CM), label, ha="center", va="top",
                 fontsize=FS_GROUP, fontweight="bold", color="0.15")
        fig.add_artist(Line2D([x0 + pad, x1 - pad], [y_of(HEADER_GROUP_RULE_CM)] * 2, transform=fig.transFigure,
                              lw=0.7, color="0.25", solid_capstyle="butt"))
    xc_all = 0.5 * (hax[0].get_position().x0 + hax[-1].get_position().x1)
    legend_h = [Patch(color="0.25", label="+ strand (dark shade)"), Patch(color="0.75", label="\u2212 strand (light shade)"),
                Line2D([], [], color="black", lw=1.4, label="+ and \u2212 pooled (grey band: 95% CI)"),
                Line2D([], [], color="red", lw=1.0, ls="--", label="ratio = 1")]
    fig.legend(handles=legend_h, loc="upper center", bbox_to_anchor=(xc_all, y_top + TOP_LEGEND_DROP),
               ncol=len(legend_h), frameon=False, fontsize=FS_LEGEND, handlelength=1.6, columnspacing=1.4,
               borderpad=0, handletextpad=0.5)

    abs_fmt = FuncFormatter(lambda v, _: f"{abs(v):g}")

    def add_guides(ax):
        """faint vertical guide lines at the major ticks, on top of the bars (reading aid)"""
        x0, x1 = sorted(ax.get_xlim())
        for t in ax.get_xticks():
            if x0 + 1e-9 < t < x1 - 1e-9 and abs(t) > 1e-9 and abs(t - 1.0) > 1e-9:
                ax.axvline(t, color="0.45", alpha=0.35, lw=0.45, zorder=5)

    for r, f in enumerate(FEATURES):
        ax0 = None                      # every row has its OWN length axis (not shared between classes)
        d = res[f]
        cdark, clight = COLORS[f]
        axes = []
        for j in range(n_cols):
            ax = fig.add_subplot(gs[r + 1, j], sharey=ax0)
            if j == 0:
                ax0 = ax
            axes.append(ax)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            ax.tick_params(axis="both", labelsize=FS_TICK)
            if r == n_rows - 1:
                ax.set_xlabel(xlabels[j], fontsize=FS_LABEL, labelpad=1.5)
            if j > 0:
                ax.tick_params(axis="y", left=False, labelleft=False)
        axes[0].annotate(f, xy=(0, 0.5), xycoords="axes fraction", xytext=(-46, 0), textcoords="offset points",
                         rotation=90, ha="center", va="center", fontsize=FS_ROW, fontweight="bold")
        if d["status"] != "ok":
            for ax in axes:
                ax.text(0.5, 0.5, "insufficient", ha="center", va="center", transform=ax.transAxes, fontsize=FS_STAT)
                ax.set_xticks([])
            axes[0].set_yticks([])
            continue
        A, W = d["arr"], d["wide"]
        y, bw, e = A["ycent"], A["bin_w"], A["edges"]

        # ---- column 1: length pattern -------------------------------------------------------------
        ax = axes[0]
        ax.stairs(A["len_p"], e, orientation="horizontal", fill=True, color=cdark, lw=0, baseline=0)
        ax.stairs(-A["len_m"], e, orientation="horizontal", fill=True, color=clight, lw=0, baseline=0)
        ax.axvline(0, color="black", lw=0.5, zorder=4)
        xm = max(A["len_p"].max(), A["len_m"].max()) * 1.08
        ax.set_xlim(-xm, xm)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=2, symmetric=True))
        ax.xaxis.set_major_formatter(abs_fmt)
        add_guides(ax)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.set_ylabel("Length (nt)", fontsize=FS_LABEL, labelpad=1.5)
        ax.set_title(label_length(W, d), fontsize=FS_STAT, pad=2.5, linespacing=1.1)

        # ---- columns 2-5: nucleotide patterns (portion of the base in every length bin) ---------------
        fam_max = {"AT": 0.0, "GC": 0.0}
        for a_, c_, *_ in PAIRS:
            fam = "AT" if a_ in (0, 1) else "GC"
            fam_max[fam] = max(fam_max[fam], float(np.nanmax(A["Pp"][:, a_])), float(np.nanmax(A["Pm"][:, c_])))
        for k, (a_, c_, key, _, _) in enumerate(PAIRS):
            ax = axes[1 + k]
            xm = fam_max["AT" if a_ in (0, 1) else "GC"] * 1.05
            ax.barh(y, -A["Pm"][:, c_], height=bw, color=clight, lw=0)
            ax.barh(y, A["Pp"][:, a_], height=bw, color=cdark, lw=0)
            ax.axvline(0, color="black", lw=0.5, zorder=4)
            ax.set_xlim(-xm, xm)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=2, symmetric=True))
            ax.xaxis.set_major_formatter(abs_fmt)
            add_guides(ax)
            ax.set_title(label_pair(W, key), fontsize=FS_STAT, pad=2.5, linespacing=1.1)

        # ---- columns 6-7: A/T and G/C ratios (+ only, - only, pooled) ---------------------------------
        for k, key in enumerate(("AT", "GC")):
            ax = axes[5 + k]
            Rp, Rm, Rt, se = A["ratio"][key]
            ax.fill_betweenx(y, Rt - CI_Z * se, Rt + CI_Z * se, color="black", alpha=CI_ALPHA, lw=0)
            ax.plot(Rp, y, color=cdark, lw=0.8)
            ax.plot(Rm, y, color=clight, lw=0.8)
            ax.plot(Rt, y, color="black", lw=0.9)
            for vals_, col_ in ((Rp, cdark), (Rm, clight), (Rt, "black")):       # dots on the lines
                ax.scatter(vals_, y, s=2.2, color=col_, lw=0, zorder=3)
            ax.axvline(1.0, color="red", ls="--", lw=0.8, alpha=0.9, zorder=4)
            vals = np.concatenate([np.abs(v[np.isfinite(v)] - 1) for v in (Rp, Rm, Rt)])
            lim = 0.1 if vals.size == 0 else round(float(np.clip(np.ceil(np.percentile(vals, 95) / 0.1) * 0.1, 0.1, 1.0)), 6)
            ax.set_xlim(1 - lim, 1 + lim)
            ax.xaxis.set_major_locator(FixedLocator([1 - lim, 1.0, 1 + lim]))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.2f}".rstrip("0").rstrip(".")))
            ax.set_title(label_ratio(W, key), fontsize=FS_STAT, pad=2.5, linespacing=1.1)
        axes[0].set_ylim(*DISPLAY_Y_RANGE.get(f, (e[0], e[-1])))
    if footnote:
        fig.text(0.5, 0.004, footnote, ha="center", va="bottom", fontsize=6.5, style="italic", color="0.35")
    out_dir = Path(out_dir)
    saved = []
    for ext, kw in (("pdf", {}), ("svg", {}), ("tiff", dict(dpi=600, pil_kwargs={"compression": "tiff_lzw"})),
                    ("png", dict(dpi=300))):
        p = out_dir / f"{base}.{ext}"
        fig.savefig(p, facecolor="white", **kw)
        saved.append(p)
    plt.close(fig)
    return saved


def matrix_grid(mats, title, out_dir, name):
    """Two rows (correlation r, distance Delta) x seven classes of 8 x 8 pattern matrices.
    Head-to-head order: the first row is at the BOTTOM and the first column at the LEFT, so A+ meets A+ at the
    bottom-left corner and the diagonal (identical patterns) runs from bottom-left to top-right."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    apply_style()
    plt.rcParams.update({"font.size": 6})
    fig, axes = plt.subplots(2, len(FEATURES), figsize=(24 / 2.54, 9.5 / 2.54), facecolor="white")
    fig.subplots_adjust(left=0.05, right=0.995, top=0.88, bottom=0.09, wspace=0.28, hspace=0.35)
    labs = [pretty(l) for l in DIST_LABELS]
    dmax = 1.0
    for f in FEATURES:
        if mats[f] is not None:
            dmax = max(dmax, float(np.nanpercentile(mats[f][1], 95)))
    for j, f in enumerate(FEATURES):
        for i in range(2):
            ax = axes[i, j]
            if mats[f] is None:
                ax.axis("off")
                continue
            M = mats[f][i]
            ax.imshow(M, origin="lower", cmap="viridis" if i == 0 else "viridis_r", vmin=-1 if i == 0 else 0,
                      vmax=1 if i == 0 else dmax, aspect="equal", interpolation="nearest")
            ax.set_xticks(range(8))
            ax.set_xticklabels(labs, fontsize=5, rotation=90)
            ax.set_yticks(range(8))
            ax.set_yticklabels(labs, fontsize=5)
            ax.tick_params(length=1.5, pad=1)
            for r_ in range(8):
                for c_ in range(8):
                    v = M[r_, c_]
                    if not np.isfinite(v):
                        continue
                    frac = (v + 1) / 2 if i == 0 else 1 - min(v / dmax, 1)
                    ax.text(c_, r_, f"{v:.2f}" if i == 0 else f"{v:.1f}", ha="center", va="center", fontsize=3.6,
                            color="white" if frac < 0.55 else "black")
            for r_ in range(8):
                c_ = partner_index(r_)
                ax.add_patch(Rectangle((c_ - 0.5, r_ - 0.5), 1, 1, fill=False, ec="#d62728", lw=0.8))
            if i == 0:
                ax.set_title(f, fontsize=8, fontweight="bold", pad=3)
    axes[0, 0].set_ylabel("pattern correlation r\n(1 = same pattern)", fontsize=6.5)
    axes[1, 0].set_ylabel("distance \u0394 (pp)\n(0 = same pattern)", fontsize=6.5)
    fig.suptitle(title, fontsize=8, y=0.985)
    paths = []
    for ext, kw in (("pdf", {}), ("png", dict(dpi=300))):
        p = Path(out_dir) / f"{name}.{ext}"
        fig.savefig(p, facecolor="white", **kw)
        paths.append(p)
    plt.close(fig)
    return paths


def draw_pattern_matrices(res, chrom, out_dir):
    """Supplementary figure for the main chromosome: 8 x 8 correlation and distance matrices per class."""
    mats = {f: ((res[f]["arr"]["R"], res[f]["arr"]["D"]) if res[f]["arr"] is not None else None) for f in FEATURES}
    return matrix_grid(mats, f"{chrom}: all eight nucleotide patterns against each other (A+ with A+ = maximum, "
                             "bottom-left to top-right); red box = RC partner of each row", out_dir,
                       "Supp_Fig_chr1_pattern_matrix")


# =============================================================================
# 9. CHROMOSOME REPORT (main mode)
# =============================================================================

def short_pair(p):
    return p.replace("(+)", "+").replace("(\u2212)", "-").replace(" vs ", " v ")


def build_chr_report(res, chrom, params_line):
    pairs = pd.DataFrame([r for f in FEATURES for r in res[f]["pairs"]])
    ratios = pd.DataFrame([r for f in FEATURES for r in res[f]["ratios"]])
    cross = pd.DataFrame([r for f in FEATURES for r in res[f]["cross"]])
    wide = pd.DataFrame([res[f]["wide"] for f in FEATURES])
    mats = []
    for f in FEATURES:
        A = res[f]["arr"]
        if A is None:
            continue
        for name, M in (("r", A["R"]), ("Delta_pp", A["D"])):
            for i in range(8):
                for j in range(8):
                    mats.append({"Feature": f, "Measure": name, "Row": DIST_LABELS[i], "Col": DIST_LABELS[j],
                                 "Value": M[i, j]})
    mats = pd.DataFrame(mats)
    ok = wide[wide["Status"] == "ok"]
    L = [f"Figure 4 statistics, {chrom} (Figure4_final v{VERSION}, {time.strftime('%Y-%m-%d %H:%M')})",
         params_line, "",
         "Patterns run along the 10-nt length bins. r = Spearman correlation of two patterns (1 = same pattern), "
         "Delta = base-weighted mean |difference| of two patterns (pp).", "",
         "1. LENGTH PATTERN  (% of blocks per length bin, + vs -; D = KS distance of the length distributions)"]
    t = pd.DataFrame({"Class": ok["Feature"], "N+": ok["N_plus"], "N-": ok["N_minus"],
                      "bins": ok["Bins_used"], "balance %": ok["Len_balance_pct"].round(1), "D": ok["Len_KS_D"].round(4),
                      "chance": ok["Len_KS_D_expected"].round(4), "x chance": ok["Len_D_over_expected"].round(2),
                      "r (bin pattern)": ok["Len_r"].round(3)})
    L += [t.to_string(index=False), "",
          "2. NUCLEOTIDE PATTERNS, the four target pairs (self = same pattern from odd vs even blocks = the ceiling for this sampling)"]
    t = pd.DataFrame({"Class": pairs["Feature"], "Pair": pairs["Pair"].map(short_pair),
                      "mean+ %": pairs["Mean_plus_pct"].round(2), "mean- %": pairs["Mean_minus_pct"].round(2),
                      "r": pairs["r"].round(3), "r self": pairs["r_self"].round(3),
                      "Delta pp": pairs["Delta_pp"].round(3), "Delta self": pairs["Delta_self_pp"].round(3),
                      "rank r (of 7)": pairs["Rank_r_among_7"], "rank Delta (of 7)": pairs["Rank_Delta_among_7"]})
    L += [t.to_string(index=False), "",
          "3. CROSS TEST  all eight patterns against all eight (diagonal r = 1, Delta = 0); rank 1 = the RC partner is the nearest"]
    for f in FEATURES:
        A = res[f]["arr"]
        if A is None:
            continue
        for name, M, nd in (("r", A["R"], 2), ("Delta (pp)", A["D"], 2)):
            t = pd.DataFrame(M.round(nd), index=DIST_LABELS, columns=DIST_LABELS)
            self_ = A["r_self"] if name == "r" else A["d_self"]
            t.insert(0, "self split", np.round(self_, 3))
            L += [f"  {f}  {name}", t.to_string(), ""]
    L.append("4. NEAREST PATTERN  rank of the RC partner among the other seven (and among the four of the other strand)")
    t = pd.DataFrame({"Class": cross["Feature"], "Row": cross["Row"], "RC partner": cross["RC_partner"],
                      "r partner": cross["r_partner"].round(3), "r best non-partner": cross["r_best_nonpartner"].round(3),
                      "Delta partner": cross["Delta_partner"].round(3),
                      "Delta best non-partner": cross["Delta_best_nonpartner"].round(3),
                      "rank r /7": cross["Rank_r_among_7"], "rank Delta /7": cross["Rank_Delta_among_7"],
                      "rank r /4": cross["Rank_r_among_4_other_strand"],
                      "rank Delta /4": cross["Rank_Delta_among_4_other_strand"]})
    L += [t.to_string(index=False), ""]
    t = pd.DataFrame({"Class": ok["Feature"], "r nearest of 7 (n/8)": ok["Cross_r_nearest_of7"],
                      "Delta nearest of 7 (n/8)": ok["Cross_d_nearest_of7"],
                      "r nearest of 4 (n/8)": ok["Cross_r_nearest_of4"],
                      "Delta nearest of 4 (n/8)": ok["Cross_d_nearest_of4"]})
    L += ["Summary per class", t.to_string(index=False), "",
          "5. BASE RATIOS per length bin (reference counts): deviation = mean |R-1|; silencing = mean(+ only, - only) / pooled",
          "   r opposite = correlation of (R+ - 1) with (R- - 1) across bins (expected -1 where the strands deviate)"]
    t = pd.DataFrame({"Class": ratios["Feature"], "Ratio": ratios["Ratio"],
                      "dev + only": ratios["Dev_plus"].round(3), "dev - only": ratios["Dev_minus"].round(3),
                      "dev pooled": ratios["Dev_pooled"].round(3), "pooled noise": ratios["Pooled_noise"].round(3),
                      "silencing x": ratios["Silencing"].round(1), "r opposite": ratios["r_opposite"].round(2),
                      "R+ median": ratios["Rplus_median"].round(3), "R- median": ratios["Rminus_median"].round(3)})
    L += [t.to_string(index=False), "",
          "How to read: where A and T (or C and G) behave alike along the length bins, several patterns are about equally near and the",
          "ranking is not informative; it is informative for classes whose strands differ in composition. Patterns of different",
          "bases share the trend with length, so also compare r with the self baseline and look at Delta (level + shape)."]
    return "\n".join(L) + "\n", pairs, ratios, wide, mats

# =============================================================================
# 10. POPULATION: per-organism worker, resume, run
# =============================================================================

def settings_hash():
    d = {"v": VERSION, "bin": BIN_WIDTH_NT, "win": {f: class_window(f) for f in FEATURES},
         "min_el": MIN_ELEMENTS_PER_STRAND_BIN, "trim": TRIM_LENGTH_RANGE, "gap": MAX_GAP_BINS,
         "maxlen": MAX_FEATURE_LENGTH, "minblocks": MIN_BLOCKS_PER_STRAND, "minchr": MIN_CHROMOSOME_BP}
    return hashlib.sha1(json.dumps(d, sort_keys=True).encode()).hexdigest()[:8]


def org_dir():
    return Path(OUTPUT_DIR) / f"per_organism_{settings_hash()}"


def _apply_settings(s):
    global MIN_CHROMOSOME_BP, CHUNK_MB, SAME_RANGE_FOR_ALL, TRIM_LENGTH_RANGE
    MIN_CHROMOSOME_BP, CHUNK_MB = s["min_chromosome_bp"], s["chunk_mb"]
    SAME_RANGE_FOR_ALL, TRIM_LENGTH_RANGE = s["same_range"], s["trim"]


def process_organism(gdir):
    org = os.path.basename(gdir.rstrip("/\\"))
    cov = {"Organism": org, "fasta_files": 0, "records_total": 0, "records_rejected_name": 0,
           "records_too_short": 0, "duplicate_ids_skipped": 0, "chromosomes_analysed": 0,
           "chromosomes_failed": 0, "bp_analysed": 0}
    cdict = defaultdict(int)
    rows, fails, seen = [], [], set()
    files = find_fasta_files(gdir)
    cov["fasta_files"] = len(files)
    print(f"[ORG] {org}  ({len(files)} FASTA file(s))", flush=True)
    for fpath in files:
        try:
            for seq_id, _h, blocks, length in iter_chromosomes(fpath, None, cdict, seen):
                t0 = time.time()
                try:
                    chrom_rows = []
                    for f in FEATURES:
                        r = analyse_feature(f, blocks[f][0], blocks[f][1], keep=False)["wide"]
                        chrom_rows.append({"Organism": org, "Chromosome": seq_id, "Length_bp": length, **r})
                    rows.extend(chrom_rows)
                    cov["chromosomes_analysed"] += 1
                    cov["bp_analysed"] += int(length)
                    print(f"  [OK] {org} {seq_id} {length / 1e6:7.1f} Mb  {time.time() - t0:5.1f} s", flush=True)
                except Exception as e:
                    cov["chromosomes_failed"] += 1
                    fails.append({"Organism": org, "Chromosome": seq_id, "Error": repr(e)[:200]})
                    print(f"  [FAIL] {org} {seq_id}: {e!r}", flush=True)
                finally:
                    blocks = None
                    gc.collect()
        except Exception as e:
            cov["chromosomes_failed"] += 1
            fails.append({"Organism": org, "Chromosome": os.path.basename(fpath), "Error": "file " + repr(e)[:200]})
            print(f"  [FAIL] {fpath}: {e!r} - rest of this file skipped", flush=True)
    for k in ("records_total", "records_rejected_name", "records_too_short", "duplicate_ids_skipped"):
        cov[k] += cdict[k]
    return org, rows, fails, cov


def run_population(workers, fresh, organisms):
    d = org_dir()
    d.mkdir(parents=True, exist_ok=True)
    dirs = get_organism_dirs(INPUT_ROOT, organisms)
    if not dirs:
        print("[ERROR] no organism folder under", INPUT_ROOT)
        return
    todo = []
    for g in dirs:
        o = safe_name(os.path.basename(g.rstrip("/\\")))
        if fresh:
            for p in (d / f"{o}.csv", d / f"{o}.cov.json"):
                if p.exists():
                    p.unlink()
        if not (d / f"{o}.cov.json").exists():
            todo.append(g)
    print(f"[INFO] {len(dirs)} organism folders, {len(dirs) - len(todo)} already done, {len(todo)} to run")
    if not todo:
        return
    settings = {"min_chromosome_bp": MIN_CHROMOSOME_BP, "chunk_mb": CHUNK_MB,
                "same_range": SAME_RANGE_FOR_ALL, "trim": TRIM_LENGTH_RANGE}
    t0 = time.time()

    def handle(res):
        org, rows, fails, cov = res
        o = safe_name(org)
        if rows:
            pd.DataFrame(rows).to_csv(d / f"{o}.csv", index=False)
        (d / f"{o}.cov.json").write_text(json.dumps(cov))
        if fails:
            p = Path(OUTPUT_DIR) / "failed_chromosomes.csv"
            pd.DataFrame(fails).to_csv(p, mode="a", header=not p.exists(), index=False)
        print(f"[DONE ORG] {org}: {cov['chromosomes_analysed']} chromosomes "
              f"({(time.time() - t0) / 60:.1f} min elapsed)", flush=True)

    if workers > 1 and len(todo) > 1:
        try:
            ctx = mp.get_context("fork")
        except ValueError:
            ctx = mp.get_context()
        with ctx.Pool(min(workers, len(todo)), initializer=_apply_settings, initargs=(settings,)) as pool:
            for res in pool.imap_unordered(process_organism, todo):
                handle(res)
    else:
        for g in todo:
            handle(process_organism(g))
    print(f"[INFO] finished in {(time.time() - t0) / 60:.1f} min")


# =============================================================================
# 11. POPULATION SUMMARY (tables + report for the Results text)
# =============================================================================

def load_results():
    files = sorted(glob.glob(str(org_dir() / "*.csv")))
    if not files:
        return None, None
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df.drop_duplicates(["Organism", "Chromosome", "Feature"], keep="last")
    df["Feature"] = pd.Categorical(df["Feature"], FEATURES, ordered=True)
    df = df.sort_values(["Organism", "Chromosome", "Feature"]).reset_index(drop=True)
    covs = [json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(org_dir() / "*.cov.json")))]
    return df, pd.DataFrame(covs)


def stacked(df, suffix, absval=False):
    """The four target pairs of every chromosome stacked into one column (Feature, Pair, v)."""
    parts = []
    for key in PAIR_KEYS:
        v = df[f"{key}_{suffix}"]
        parts.append(pd.DataFrame({"Feature": df["Feature"].astype(str), "Pair": key, "v": v.abs() if absval else v}))
    return pd.concat(parts, ignore_index=True)


def stacked_ratio(df, suffix):
    parts = []
    for key in ("AT", "GC"):
        parts.append(pd.DataFrame({"Feature": df["Feature"].astype(str), "Ratio": key, "v": df[f"{key}_{suffix}"]}))
    return pd.concat(parts, ignore_index=True)


def stats5(s):
    s = pd.Series(s).dropna().astype(float)
    if s.empty:
        return [np.nan] * 5
    return [s.median(), s.quantile(.25), s.quantile(.75), s.min(), s.max()]


def fmt5(s, d=3):
    m, q1, q3, lo, hi = stats5(s)
    if not np.isfinite(m):
        return "n/a"
    return f"{m:.{d}f} (IQR {q1:.{d}f}-{q3:.{d}f}; range {lo:.{d}f}-{hi:.{d}f})"


def rank_cols(sub, kind):
    return pd.concat([sub[f"Rank{kind}_{l}"] for l in DIST_LABELS], ignore_index=True)


def pct_rank1(sub, kind):
    r = rank_cols(sub, kind).dropna()
    return 100 * float((r == 1).mean()) if len(r) else np.nan


def pct_of4(sub, col):
    return 100 * float(sub[col].sum()) / (8 * len(sub)) if len(sub) else np.nan


def class_summary(ok):
    rows = []
    for f in FEATURES:
        s = ok[ok["Feature"] == f]
        if s.empty:
            continue
        r = {"Feature": f, "Chromosomes": len(s), "Organisms": s["Organism"].nunique()}
        series = (("Len_D", s["Len_KS_D"]), ("Len_D_ratio", s["Len_D_over_expected"]), ("Len_r", s["Len_r"]),
                  ("Pairs_r", stacked(s, "r")["v"]), ("Pairs_r_self", stacked(s, "r_self")["v"]),
                  ("Pairs_Delta_pp", stacked(s, "delta")["v"]), ("Pairs_Delta_self_pp", stacked(s, "delta_self")["v"]),
                  ("Ratio_dev_plus", stacked_ratio(s, "dev_plus")["v"]), ("Ratio_dev_minus", stacked_ratio(s, "dev_minus")["v"]),
                  ("Ratio_dev_pooled", stacked_ratio(s, "dev_pooled")["v"]),
                  ("Ratio_noise", stacked_ratio(s, "noise")["v"]), ("Ratio_silencing", stacked_ratio(s, "silencing")["v"]))
        for name, ser in series:
            for lab, v in zip(("median", "q25", "q75", "min", "max"), stats5(ser)):
                r[f"{name}_{lab}"] = v
        r["Cross_pct_nearest_r_of7"] = pct_rank1(s, "R7")
        r["Cross_pct_nearest_delta_of7"] = pct_rank1(s, "D7")
        r["Cross_pct_nearest_r_of4"] = pct_of4(s, "Cross_r_nearest_of4")
        r["Cross_pct_nearest_delta_of4"] = pct_of4(s, "Cross_d_nearest_of4")
        r["Cross_pct_chr_all8_delta_of4"] = 100 * float((s["Cross_d_nearest_of4"] == 8).mean())
        for k in ("r_partner", "r_self", "r_nonpartner", "d_partner", "d_self", "d_nonpartner"):
            r[f"Cross_{k}_median"] = s[f"Cross_{k}_median"].median()
        rows.append(r)
    return pd.DataFrame(rows)


def build_population_report(df, ok, covdf):
    n_org, n_chr = df["Organism"].nunique(), df.drop_duplicates(["Organism", "Chromosome"]).shape[0]
    L = [f"Population check of + / - annotated populations (Figure4_final v{VERSION}, {time.strftime('%Y-%m-%d %H:%M')})", "",
         f"Organisms: {n_org}   Chromosomes: {n_chr}   Chromosome x class analyses: {len(df)}   "
         f"usable: {len(ok)}   insufficient (< {MIN_BLOCKS_PER_STRAND} blocks per strand or < {MIN_USABLE_BINS} usable length bins): "
         f"{len(df) - len(ok)}"]
    if covdf is not None and len(covdf):
        c = covdf.sum(numeric_only=True)
        L.append(f"FASTA records read: {int(c['records_total'])}; rejected by name: {int(c['records_rejected_name'])}; "
                 f"shorter than {MIN_CHROMOSOME_BP:,} bp: {int(c['records_too_short'])}; "
                 f"duplicate ids skipped: {int(c['duplicate_ids_skipped'])}; failed: {int(c['chromosomes_failed'])}")
        none = covdf[covdf["chromosomes_analysed"] == 0]["Organism"].tolist()
        if none:
            L.append(f"Organisms with NO analysed chromosome ({len(none)}): " + ", ".join(none[:15])
                     + (" ..." if len(none) > 15 else ""))
    ins = df[df["Status"] != "ok"].groupby("Feature", observed=True).size()
    L.append("Insufficient per feature: " + ", ".join(f"{f} {int(ins.get(f, 0))}" for f in FEATURES))
    L += ["", f"Length bins {BIN_WIDTH_NT} nt; >= {MIN_ELEMENTS_PER_STRAND_BIN} blocks per strand per bin; tail trimming "
              f"{'on' if TRIM_LENGTH_RANGE else 'off'} (gap {MAX_GAP_BINS}); no RC transformation of the - data; nothing shuffled.",
          "Unit = one chromosome x one class. Values are median (IQR; range) across chromosomes.",
          "Patterns run along the length bins. r = Spearman correlation of two patterns (1 = same pattern); Delta = base-weighted",
          "mean |difference| of two patterns (pp); D = KS distance of two histograms; x chance = D / D expected by sampling alone.", ""]
    cs = class_summary(ok)
    L.append("A. LENGTH PATTERN  (% of blocks per length bin, + vs -)")
    for _, r in cs.iterrows():
        s = ok[ok["Feature"] == r["Feature"]]
        L.append(f"  {r['Feature']:<7} chr {int(r['Chromosomes']):>5}  D {fmt5(s['Len_KS_D'], 3)}  x chance {fmt5(s['Len_D_over_expected'], 2)}")
        L.append(f"          r of the bin pattern {fmt5(s['Len_r'], 3)}   balance % {s['Len_balance_pct'].median():.1f}")
    L += ["", "B. NUCLEOTIDE PATTERNS, target pairs (A+ vs T-, T+ vs A-, G+ vs C-, C+ vs G-; four pairs pooled)"]
    for _, r in cs.iterrows():
        s = ok[ok["Feature"] == r["Feature"]]
        L.append(f"  {r['Feature']:<7} r {fmt5(stacked(s, 'r')['v'], 3)}   self-split r {stacked(s, 'r_self')['v'].median():.3f}")
        L.append(f"          Delta pp {fmt5(stacked(s, 'delta')['v'], 3)}   self-split Delta {stacked(s, 'delta_self')['v'].median():.3f}")
    L += ["", "   per pair, all classes pooled:"]
    for key, (_, _, _, _, name) in zip(PAIR_KEYS, PAIRS):
        L.append(f"  {name:<14} r {fmt5(ok[f'{key}_r'], 3)}   Delta pp {fmt5(ok[f'{key}_delta'], 3)}")
    L += ["", "C. CROSS TEST AND NEAREST PATTERN  (eight patterns per class; RC partner = T- for A+, A- for T+, ...; rank 1 = nearest)"]
    for _, r in cs.iterrows():
        L.append(f"  {r['Feature']:<7} partner nearest of 7: by r {r['Cross_pct_nearest_r_of7']:.1f}%, by Delta {r['Cross_pct_nearest_delta_of7']:.1f}%"
                 f";  of 4 (other strand): by r {r['Cross_pct_nearest_r_of4']:.1f}%, by Delta {r['Cross_pct_nearest_delta_of4']:.1f}%"
                 f";  all 8 rows right (Delta, of 4) in {r['Cross_pct_chr_all8_delta_of4']:.1f}% of chromosomes")
        L.append(f"          median r: partner {r['Cross_r_partner_median']:.3f} | self {r['Cross_r_self_median']:.3f} | best non-partner "
                 f"{r['Cross_r_nonpartner_median']:.3f};   median Delta: partner {r['Cross_d_partner_median']:.3f} | self "
                 f"{r['Cross_d_self_median']:.3f} | best non-partner {r['Cross_d_nonpartner_median']:.3f}")
    r7r, r7d = rank_cols(ok, "R7").dropna(), rank_cols(ok, "D7").dropna()
    L += ["", f"  All classes pooled: partner nearest of 7 by r in {100 * (r7r == 1).mean():.1f}% and by Delta in "
              f"{100 * (r7d == 1).mean():.1f}% of {len(r7d):,} rows.", "",
          "D. BASE RATIOS per length bin (reference counts): deviation = mean |R-1|; silencing = mean(+ only, - only) / pooled"]
    for _, r in cs.iterrows():
        s = ok[ok["Feature"] == r["Feature"]]
        L.append(f"  {r['Feature']:<7} dev + only {stacked_ratio(s, 'dev_plus')['v'].median():.3f} | - only "
                 f"{stacked_ratio(s, 'dev_minus')['v'].median():.3f} | pooled {stacked_ratio(s, 'dev_pooled')['v'].median():.3f} "
                 f"(noise {stacked_ratio(s, 'noise')['v'].median():.3f});  silencing x {fmt5(stacked_ratio(s, 'silencing')['v'], 1)}")
    L += ["", "FOR THE RESULTS TEXT (all chromosomes, all classes pooled):",
          f"  length D                     {fmt5(ok['Len_KS_D'], 3)}",
          f"  length x chance              {fmt5(ok['Len_D_over_expected'], 2)}",
          f"  length bin-pattern r         {fmt5(ok['Len_r'], 3)}",
          f"  nucleotide pattern r         {fmt5(stacked(ok, 'r')['v'], 3)}",
          f"  nucleotide pattern Delta pp  {fmt5(stacked(ok, 'delta')['v'], 3)}",
          f"  RC partner nearest (of 7)    by r {100 * (r7r == 1).mean():.1f}%, by Delta {100 * (r7d == 1).mean():.1f}%",
          f"  ratio deviation + only       {fmt5(stacked_ratio(ok, 'dev_plus')['v'], 3)}",
          f"  ratio deviation - only       {fmt5(stacked_ratio(ok, 'dev_minus')['v'], 3)}",
          f"  ratio deviation pooled       {fmt5(stacked_ratio(ok, 'dev_pooled')['v'], 3)}",
          f"  silencing factor             {fmt5(stacked_ratio(ok, 'silencing')['v'], 1)}", "",
          "How to read: where A and T (or C and G) behave alike along the length bins, several patterns are about equally near and the",
          "ranking is not informative. All patterns share a trend with length, so compare r with the self-split baseline and read Delta",
          "(level + shape). Length-pattern r is high almost by construction because both strands share the same decay with length;",
          "D and x chance are the stricter length statistics. With 10^4-10^5 blocks tiny differences are detectable: read the sizes."]
    fp = Path(OUTPUT_DIR) / "failed_chromosomes.csv"
    if fp.exists():
        L += ["", f"WARNING: some chromosomes/files failed, see {fp.name}"]
    return "\n".join(L) + "\n", cs


# =============================================================================
# 12. SUPPLEMENTARY FIGURES (population)
# =============================================================================

def _setup_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    apply_style()
    plt.rcParams.update({"font.size": 7})
    return plt


def _clean(a):
    a = np.asarray(a, float)
    return a[np.isfinite(a)]


def _boxes(ax, data, pos, color, width=0.7):
    keep = [(_clean(d), p) for d, p in zip(data, pos) if len(_clean(d))]
    if not keep:
        return
    bp = ax.boxplot([k[0] for k in keep], positions=[k[1] for k in keep], widths=width, patch_artist=True,
                    showfliers=True, flierprops=dict(marker=".", markersize=1.5, markerfacecolor=color,
                                                     markeredgecolor=color, alpha=0.35),
                    medianprops=dict(color="black", lw=0.8), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5),
                    boxprops=dict(lw=0.5))
    for b in bp["boxes"]:
        b.set_facecolor(color)
        b.set_alpha(0.85)


def _style(ax, letter, title, ylabel=None):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(title, fontsize=7.5, fontweight="bold", pad=4)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=7)
    ax.text(-0.14, 1.12, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


def _feat_ticks(ax, step=1, offset=0.0):
    ax.set_xticks(np.arange(len(FEATURES)) * step + offset)
    ax.set_xticklabels(FEATURES, rotation=40, ha="right", fontsize=6.5)


def _save(fig, out_dir, name):
    paths = []
    for ext, kw in (("pdf", {}), ("png", dict(dpi=300))):
        p = Path(out_dir) / f"{name}.{ext}"
        fig.savefig(p, facecolor="white", **kw)
        paths.append(p)
    return paths


def _scatter_by_class(ax, ok, xcol, ycol):
    for f in FEATURES:
        s = ok[ok["Feature"] == f]
        ax.scatter(s[xcol], s[ycol], s=3, color=COLORS[f][0], alpha=0.5, lw=0, label=f)


def _identity(ax, lo, hi):
    ax.plot([lo, hi], [lo, hi], color="black", lw=0.6)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)


def _title_suffix(ok):
    return (f"{ok.drop_duplicates(['Organism', 'Chromosome']).shape[0]} chromosomes of {ok['Organism'].nunique()} organisms "
            "(each point or box value = one chromosome of one class)")


def supp_s1_length(ok, out_dir):
    plt = _setup_axes()
    from matplotlib.ticker import FuncFormatter
    fig, axs = plt.subplots(2, 2, figsize=(17.5 / 2.54, 14 / 2.54), facecolor="white")
    fig.subplots_adjust(left=0.09, right=0.99, top=0.93, bottom=0.14, wspace=0.28, hspace=0.62)
    sub = {f: ok[ok["Feature"] == f] for f in FEATURES}
    a = axs[0, 0]
    _scatter_by_class(a, ok, "N_plus", "N_minus")
    a.set_xscale("log")
    a.set_yscale("log")
    lim = [ok[["N_plus", "N_minus"]].min().min() / 1.3, ok[["N_plus", "N_minus"]].max().max() * 1.3]
    _identity(a, *lim)
    a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    a.set_xlabel("blocks on + strand (N+)")
    _style(a, "A", "Size of the two populations", "blocks on \u2212 strand (N\u2212)")
    a.legend(fontsize=5.5, frameon=False, markerscale=2.5, handletextpad=0.2, borderpad=0)
    for a, col, letter, title, yl, hl in ((axs[0, 1], "Len_KS_D", "B", "Length distribution: + vs \u2212", "KS distance D (0 = identical)", 0),
                                          (axs[1, 0], "Len_D_over_expected", "C", "D relative to chance", "D / D expected (1 = chance)", 1),
                                          (axs[1, 1], "Len_r", "D", "Correlation of the bin patterns", "Spearman r of % blocks per length bin", None)):
        for i, f in enumerate(FEATURES):
            _boxes(a, [sub[f][col]], [i], COLORS[f][0])
        _feat_ticks(a)
        if hl is not None:
            a.axhline(hl, color="0.4", lw=0.6, ls="--" if hl == 1 else "-")
        _style(a, letter, title, yl)
    fig.suptitle("Length-pattern similarity of + and \u2212 blocks: " + _title_suffix(ok), fontsize=8, y=0.995)
    paths = _save(fig, out_dir, "Supp_Fig_S3_length_detail")
    plt.close(fig)
    return paths


def supp_s2_pairs(ok, out_dir):
    plt = _setup_axes()
    from matplotlib.patches import Patch
    fig, axs = plt.subplots(2, 3, figsize=(24 / 2.54, 14.5 / 2.54), facecolor="white")
    fig.subplots_adjust(left=0.06, right=0.99, top=0.93, bottom=0.12, wspace=0.3, hspace=0.7)
    sub = {f: ok[ok["Feature"] == f] for f in FEATURES}

    def grouped(a, suffix, letter, title, ylabel, hline=None):
        for i, f in enumerate(FEATURES):
            for j, key in enumerate(PAIR_KEYS):
                _boxes(a, [sub[f][f"{key}_{suffix}"]], [i * 5 + j], PAIR_COLORS[key], width=0.8)
        _feat_ticks(a, 5, 1.5)
        if hline is not None:
            a.axhline(hline, color="0.5", lw=0.5)
        _style(a, letter, title, ylabel)

    grouped(axs[0, 0], "r", "A", "Nucleotide pattern: correlation of the RC pair", "Spearman r across length bins", hline=0)
    grouped(axs[0, 1], "delta", "B", "Nucleotide pattern: distance of the RC pair", "\u0394 (pp)", hline=0)
    a = axs[0, 2]
    for key in PAIR_KEYS:
        a.scatter(ok[f"{key}_r_self"], ok[f"{key}_r"], s=3, color=PAIR_COLORS[key], alpha=0.4, lw=0)
    _identity(a, -0.2, 1.02)
    a.set_xlabel("r of the same pattern, odd vs even blocks (self)")
    _style(a, "C", "RC pair against the self baseline: r", "r of the RC pair")
    a.text(0.03, 0.95, "on the line = as similar as itself", transform=a.transAxes, ha="left", va="top", fontsize=6)
    a = axs[1, 0]
    for key in PAIR_KEYS:
        a.scatter(ok[f"{key}_delta_self"], ok[f"{key}_delta"], s=3, color=PAIR_COLORS[key], alpha=0.4, lw=0)
    mx = float(np.nanpercentile(np.concatenate([ok[f"{k}_delta"].dropna().to_numpy() for k in PAIR_KEYS]), 99.5)) * 1.05
    _identity(a, 0, mx)
    a.set_xlabel("\u0394 of the same pattern, odd vs even blocks (self)")
    _style(a, "D", "RC pair against the self baseline: \u0394", "\u0394 of the RC pair (pp)")
    a = axs[1, 1]
    for key in PAIR_KEYS:
        a.scatter(ok[f"{key}_mean_minus"], ok[f"{key}_mean_plus"], s=3, color=PAIR_COLORS[key], alpha=0.4, lw=0)
    allm = np.concatenate([ok[f"{k}_mean_plus"].dropna().to_numpy() for k in PAIR_KEYS]
                          + [ok[f"{k}_mean_minus"].dropna().to_numpy() for k in PAIR_KEYS])
    _identity(a, float(allm.min()) - 1, float(allm.max()) + 1)
    a.set_xlabel("portion of the paired base in \u2212 blocks (%)")
    _style(a, "E", "Mean portion: + against its RC partner in \u2212", "portion of the base in + blocks (%)")
    axs[1, 2].axis("off")
    axs[1, 2].legend(handles=[Patch(color=PAIR_COLORS[k], label=p[4]) for k, p in zip(PAIR_KEYS, PAIRS)],
                     loc="center", frameon=False, fontsize=8, title="target pairs", title_fontsize=8)
    fig.suptitle("Nucleotide-pattern similarity of + and \u2212 blocks (target pairs): " + _title_suffix(ok), fontsize=8, y=0.995)
    paths = _save(fig, out_dir, "Supp_Fig_S4_pairs_detail")
    plt.close(fig)
    return paths


def supp_s3_cross(ok, out_dir):
    plt = _setup_axes()
    fig, axs = plt.subplots(2, 2, figsize=(17.5 / 2.54, 14 / 2.54), facecolor="white")
    fig.subplots_adjust(left=0.09, right=0.99, top=0.93, bottom=0.14, wspace=0.28, hspace=0.62)
    sub = {f: ok[ok["Feature"] == f] for f in FEATURES}
    x = np.arange(len(FEATURES))
    a = axs[0, 0]
    bars = (("by r, of 7", [pct_rank1(sub[f], "R7") for f in FEATURES], "#0072B2"),
            ("by \u0394, of 7", [pct_rank1(sub[f], "D7") for f in FEATURES], "#56B4E9"),
            ("by r, of 4 (other strand)", [pct_of4(sub[f], "Cross_r_nearest_of4") for f in FEATURES], "#E69F00"),
            ("by \u0394, of 4 (other strand)", [pct_of4(sub[f], "Cross_d_nearest_of4") for f in FEATURES], "#D55E00"))
    for k, (lab, vals, col) in enumerate(bars):
        a.bar(x - 0.3 + 0.2 * k, vals, 0.19, color=col, label=lab)
    a.set_ylim(0, 112)
    a.set_xticks(x)
    a.set_xticklabels(FEATURES, rotation=40, ha="right", fontsize=6.5)
    _style(a, "A", "RC partner is the nearest pattern", "% of patterns (class x chromosome x 8)")
    a.legend(fontsize=5.2, frameon=False, loc="upper center", ncol=2, columnspacing=0.8)
    a = axs[0, 1]
    bottom = np.zeros(len(FEATURES))
    cols = {"rank 1": "#009E73", "rank 2": "#E69F00", "rank 3\u20137": "#D55E00"}
    shares = {k: [] for k in cols}
    for f in FEATURES:
        r7 = rank_cols(sub[f], "D7").dropna() if len(sub[f]) else pd.Series(dtype=float)
        n = max(len(r7), 1)
        shares["rank 1"].append(100 * float((r7 == 1).sum()) / n)
        shares["rank 2"].append(100 * float((r7 == 2).sum()) / n)
        shares["rank 3\u20137"].append(100 * float((r7 >= 3).sum()) / n)
    for k, c in cols.items():
        a.bar(x, shares[k], 0.6, bottom=bottom, color=c, label=k)
        bottom += np.array(shares[k])
    a.set_ylim(0, 100)
    a.set_xticks(x)
    a.set_xticklabels(FEATURES, rotation=40, ha="right", fontsize=6.5)
    _style(a, "B", "Rank of the RC partner (by \u0394, among 7)", "% of patterns")
    a.legend(fontsize=5.5, frameon=False, loc="lower left", ncol=3)
    a = axs[1, 0]
    _scatter_by_class(a, ok, "Cross_r_nonpartner_median", "Cross_r_partner_median")
    _identity(a, -0.5, 1.02)
    a.set_xlabel("r of the best non-partner pattern")
    _style(a, "C", "RC partner against the best other pattern: r", "r of the RC partner")
    a.text(0.03, 0.95, "above the line = partner more similar", transform=a.transAxes, ha="left", va="top", fontsize=6)
    a = axs[1, 1]
    _scatter_by_class(a, ok, "Cross_d_partner_median", "Cross_d_nonpartner_median")
    a.set_xscale("log")
    a.set_yscale("log")
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    for axis in (a.xaxis, a.yaxis):
        axis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0), numticks=8))
        axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        axis.set_minor_formatter(NullFormatter())
    lo = max(1e-3, float(ok["Cross_d_partner_median"].min()) / 1.5)
    hi = max(float(ok["Cross_d_nonpartner_median"].max()), float(ok["Cross_d_partner_median"].max())) * 1.3
    _identity(a, lo, hi)
    a.set_xlabel("\u0394 to the RC partner (pp)")
    _style(a, "D", "RC partner against the nearest other pattern: \u0394", "\u0394 to the nearest non-partner (pp)")
    a.text(0.03, 0.95, "above the line = partner nearer", transform=a.transAxes, ha="left", va="top", fontsize=6)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], marker="o", ls="", color=COLORS[f][0], ms=4, label=f) for f in FEATURES],
               loc="lower center", ncol=len(FEATURES), frameon=False, fontsize=6.5, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Cross test and nearest pattern: " + _title_suffix(ok), fontsize=8, y=0.995)
    paths = _save(fig, out_dir, "Supp_Fig_S5_cross_nearest_detail")
    plt.close(fig)
    return paths


def supp_s4_ratios(ok, out_dir):
    plt = _setup_axes()
    from matplotlib.patches import Patch
    fig, axs = plt.subplots(2, 2, figsize=(17.5 / 2.54, 14.5 / 2.54), facecolor="white")
    fig.subplots_adjust(left=0.09, right=0.99, top=0.93, bottom=0.16, wspace=0.28, hspace=0.62)
    sub = {f: ok[ok["Feature"] == f] for f in FEATURES}
    cols = {"plus": "#0072B2", "minus": "#56B4E9", "pooled": "#000000"}
    for a, key, letter in ((axs[0, 0], "AT", "A"), (axs[0, 1], "GC", "B")):
        name = "A/T" if key == "AT" else "G/C"
        for i, f in enumerate(FEATURES):
            for j, kind in enumerate(("plus", "minus", "pooled")):
                _boxes(a, [sub[f][f"{key}_dev_{kind}"]], [i * 4 + j], cols[kind], width=0.8)
        _feat_ticks(a, 4, 1.0)
        a.axhline(0, color="0.5", lw=0.5)
        _style(a, letter, f"{name} ratio deviation per length bin", "mean |R \u2212 1|")
    a = axs[1, 0]
    for i, f in enumerate(FEATURES):
        _boxes(a, [sub[f]["AT_silencing"]], [i * 3], "#CC79A7", width=0.8)
        _boxes(a, [sub[f]["GC_silencing"]], [i * 3 + 1], "#009E73", width=0.8)
    _feat_ticks(a, 3, 0.5)
    a.set_yscale("log")
    a.axhline(1, color="0.4", lw=0.6, ls="--")
    _style(a, "C", "Silencing when + and \u2212 are combined", "mean(+ only, \u2212 only) / pooled deviation")
    a = axs[1, 1]
    for key, col in (("AT", "#CC79A7"), ("GC", "#009E73")):
        a.scatter(ok[f"{key}_noise"], ok[f"{key}_dev_pooled"], s=3, color=col, alpha=0.4, lw=0)
    mx = float(np.nanpercentile(np.r_[ok["AT_dev_pooled"].dropna(), ok["GC_dev_pooled"].dropna()], 99.5)) * 1.05
    _identity(a, 0, mx)
    a.set_xlabel("pooled deviation expected from sampling noise")
    _style(a, "D", "Pooled ratio against sampling noise", "observed pooled |R \u2212 1|")
    a.text(0.97, 0.05, "on the line = noise only", transform=a.transAxes, ha="right", fontsize=6)
    fig.legend(handles=[Patch(color=cols["plus"], label="+ only (A+/T+, G+/C+)"), Patch(color=cols["minus"], label="\u2212 only"),
                        Patch(color=cols["pooled"], label="pooled"), Patch(color="#CC79A7", label="silencing A/T"),
                        Patch(color="#009E73", label="silencing G/C")], loc="lower center", ncol=5, frameon=False,
               fontsize=6.5, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Base-ratio deviation and silencing: " + _title_suffix(ok), fontsize=8, y=0.995)
    paths = _save(fig, out_dir, "Supp_Fig_S6_ratios_detail")
    plt.close(fig)
    return paths


def supp_s1_population_7col(ok, out_dir):
    """Population version of Figure 4: the same seven columns, three statistics per column.
    One box = the chromosomes of one class (colour = class); nothing is shuffled."""
    plt = _setup_axes()
    from matplotlib.patches import Patch
    fig, axs = plt.subplots(3, 7, figsize=(24 / 2.54, 15 / 2.54), facecolor="white")
    fig.subplots_adjust(left=0.05, right=0.995, top=0.85, bottom=0.12, wspace=0.62, hspace=0.35)
    sub = {f: ok[ok["Feature"] == f] for f in FEATURES}
    n_chr = ok.drop_duplicates(["Organism", "Chromosome"]).shape[0]

    def pct_rank1_pair(s_, key):
        v = s_[f"{key}_rank_d"].dropna()
        return 100 * float((v == 1).mean()) if len(v) else np.nan

    cols = [("Length: + vs \u2212", [
        ("KS distance D", lambda s_: s_["Len_KS_D"], "box", 0, False),
        ("D / chance", lambda s_: s_["Len_D_over_expected"], "box", 1, False),
        ("r of the bin pattern", lambda s_: s_["Len_r"], "box", None, False)])]
    for key, (_, _, _, _, name) in zip(PAIR_KEYS, PAIRS):
        cols.append((name, [
            ("\u0394 of the patterns (pp)", lambda s_, k=key: s_[f"{k}_delta"], "box", 0, False),
            ("r of the patterns", lambda s_, k=key: s_[f"{k}_r"], "box", 0, False),
            ("% RC nearest", lambda s_, k=key: pct_rank1_pair(s_, k), "bar", None, False)]))
    for key, name in (("AT", "A/T ratio"), ("GC", "G/C ratio")):
        cols.append((name, [
            ("strand-only |R\u22121|", lambda s_, k=key: 0.5 * (s_[f"{k}_dev_plus"] + s_[f"{k}_dev_minus"]), "box", 0, False),
            ("pooled |R\u22121|", lambda s_, k=key: s_[f"{k}_dev_pooled"], "box", 0, False),
            ("silencing (\u00f7)", lambda s_, k=key: s_[f"{k}_silencing"], "box", 1, True)]))
    x = np.arange(len(FEATURES))
    for j, (title, rows) in enumerate(cols):
        for i, (ylabel, getter, kind, hline, logy) in enumerate(rows):
            ax = axs[i, j]
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if kind == "box":
                for k, f in enumerate(FEATURES):
                    _boxes(ax, [getter(sub[f])], [k], COLORS[f][0], width=0.7)
            else:
                ax.bar(x, [getter(sub[f]) if len(sub[f]) else np.nan for f in FEATURES], 0.7,
                       color=[COLORS[f][0] for f in FEATURES])
                ax.set_ylim(0, 105)
            if logy:
                ax.set_yscale("log")
            if hline is not None:
                ax.axhline(hline, color="0.4", lw=0.6, ls="--" if hline == 1 else "-")
            ax.set_ylabel(ylabel, fontsize=5.8, labelpad=1.5)
            ax.tick_params(axis="y", labelsize=5.5, length=2)
            ax.set_xlim(-0.7, len(FEATURES) - 0.3)
            if i == 2:
                ax.set_xticks(x)
                ax.set_xticklabels(FEATURES, rotation=90, fontsize=5.5)
                ax.tick_params(axis="x", length=2)
            else:
                ax.set_xticks([])
            if i == 0:
                ax.set_title(title, fontsize=7, fontweight="bold", pad=4)
    for i, letter in enumerate("ABC"):
        fig.text(0.003, 0.80 - i * 0.30, letter, fontsize=11, fontweight="bold")
    fig.legend(handles=[Patch(color=COLORS[f][0], label=f) for f in FEATURES], loc="upper center", ncol=len(FEATURES),
               frameon=False, fontsize=7, bbox_to_anchor=(0.5, 0.965))
    fig.suptitle(f"Population analysis of the + and \u2212 annotated populations: {n_chr} chromosomes of "
                 f"{ok['Organism'].nunique()} organisms (one box = all chromosomes of one class)", fontsize=8, y=0.995)
    paths = _save(fig, out_dir, "Supp_Fig_S1_population_7col")
    plt.close(fig)
    return paths


def supp_s2_population_matrix(ok, out_dir):
    """Median 8 x 8 pattern matrices over all chromosomes, head-to-head order (A+ with A+ bottom-left)."""
    mats = {}
    for f in FEATURES:
        s_ = ok[ok["Feature"] == f]
        if s_.empty:
            mats[f] = None
            continue
        R, D = np.eye(8), np.zeros((8, 8))
        for i in range(8):
            for j in range(i + 1, 8):
                R[i, j] = R[j, i] = s_[f"R_{DKEYS[i]}_{DKEYS[j]}"].median()
                D[i, j] = D[j, i] = s_[f"D_{DKEYS[i]}_{DKEYS[j]}"].median()
        mats[f] = (R, D)
    n_chr = ok.drop_duplicates(["Organism", "Chromosome"]).shape[0]
    return matrix_grid(mats, f"Median of {n_chr} chromosomes: all eight nucleotide patterns against each other "
                             "(A+ with A+ = maximum, bottom-left to top-right); red box = RC partner of each row",
                       out_dir, "Supp_Fig_S2_population_cross_matrix")


def write_excel(path, sheets, readme):
    """Several tables in one .xlsx (README sheet first). Skipped with a message if openpyxl is missing."""
    try:
        from openpyxl.styles import Alignment, Font
    except ImportError:
        print("[INFO] openpyxl is not installed: Excel tables skipped (the CSV files are still written)")
        return None
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        pd.DataFrame({"README": readme}).to_excel(xw, sheet_name="README", index=False)
        for name, df in sheets.items():
            if df is None or len(df) == 0:
                continue
            df = df.copy()
            for c in df.columns:
                if str(df[c].dtype) == "category":
                    df[c] = df[c].astype(str)
            df = df.round(4)
            df.to_excel(xw, sheet_name=name[:31], index=False)
        for ws in xw.book.worksheets:
            ws.freeze_panes = "A2"
            for cell in ws[1]:
                cell.font = Font(bold=True)
            if ws.title == "README":
                ws.column_dimensions["A"].width = 130
                for row in ws.iter_rows(min_row=2):
                    row[0].alignment = Alignment(wrap_text=True, vertical="top")
                continue
            for col in ws.columns:
                width = max(len(str(c.value)) if c.value is not None else 0 for c in list(col)[:200]) + 2
                ws.column_dimensions[col[0].column_letter].width = min(max(width, 8), 30)
    return Path(path)


README_COMMON = [
    "Patterns run along the 10-nt length bins; bases of all blocks of a bin are pooled; the - blocks are used in reference "
    "coordinates (no reverse complement); nothing is shuffled.",
    "D = KS distance between two histograms (0 = identical); chance / x chance = D expected for two samples of the same sizes "
    "from one distribution (0.87*sqrt((n+m)/(n*m))) and D divided by it.",
    "Delta (pp) = base-weighted mean |difference| of two nucleotide patterns over the usable length bins, in percentage points.",
    "r = Spearman correlation of two patterns across the length bins (1 = same pattern); r_self / Delta_self = the same pattern built "
    "from odd-numbered and even-numbered blocks (what an identical pattern gives with this sampling).",
    "Rank = rank of the RC partner (A+ -> T-, T+ -> A-, G+ -> C-, C+ -> G-) among the other seven patterns (1 = nearest); "
    "'worst of two views' for a pair means the larger of the ranks seen from the + and from the - pattern.",
    "Ratio deviation = base-weighted mean |R-1| per length bin for A+/T+ (+ only), A-/T- (- only) and the pooled ratio "
    "(A+ + A-)/(T+ + T-); silencing = mean(+ only, - only) / pooled; noise = pooled deviation expected from sampling "
    "(0.80 x cluster-robust standard error); r_opposite = correlation of (R+ - 1) with (R- - 1) across bins.",
]


def cross_df_for_excel(res):
    return pd.DataFrame([r for f in FEATURES for r in res[f]["cross"]])


def write_chr_excel(out_dir, tag, chrom, wide, pairs, ratios, cross, mats):
    len_cols = [c for c in ("Feature", "Status", "N_plus", "N_minus", "Bins_used", "Len_balance_pct", "Len_KS_D",
                            "Len_KS_D_expected", "Len_D_over_expected", "Len_r") if c in wide.columns]
    readme = [f"Supplementary tables for {chrom} (Figure4_final v{VERSION})", "",
              "Sheet S1a length: length distribution (N+, N-, KS D, chance, r of the bin pattern).",
              "Sheet S1b nucleotide patterns: the four target pairs (means, r, Delta, self-split baseline, ranks).",
              "Sheet S1c cross test: every one of the eight patterns against its RC partner and the best non-partner.",
              "Sheet S1d ratios: + only, - only and pooled A/T and G/C deviations, noise and silencing.",
              "Sheet S1e matrices: all 8 x 8 correlations (r) and distances (Delta, pp), long format.", ""] + README_COMMON
    return write_excel(Path(out_dir) / f"Supp_Tables_{tag}.xlsx",
                       {"S1a length": wide[len_cols], "S1b nucleotide patterns": pairs, "S1c cross test": cross,
                        "S1d ratios": ratios, "S1e matrices": mats}, readme)


def organism_summary(ok):
    t = ok.copy()
    t["Pairs_delta_mean"] = t[[f"{k}_delta" for k in PAIR_KEYS]].mean(axis=1)
    t["Pairs_r_mean"] = t[[f"{k}_r" for k in PAIR_KEYS]].mean(axis=1)
    t["RC_nearest_pct"] = 100 * t["Cross_d_nearest_of4"] / 8
    out = t.groupby(["Organism", "Feature"]).agg(
        Chromosomes=("Chromosome", "nunique"), Length_D_median=("Len_KS_D", "median"),
        Length_D_over_chance_median=("Len_D_over_expected", "median"),
        Pairs_Delta_pp_median=("Pairs_delta_mean", "median"), Pairs_r_median=("Pairs_r_mean", "median"),
        RC_nearest_pct_of4=("RC_nearest_pct", "mean"), AT_pooled_dev_median=("AT_dev_pooled", "median"),
        GC_pooled_dev_median=("GC_dev_pooled", "median"), AT_silencing_median=("AT_silencing", "median"),
        GC_silencing_median=("GC_silencing", "median")).reset_index()
    out["Feature"] = pd.Categorical(out["Feature"], FEATURES, ordered=True)
    return out.sort_values(["Organism", "Feature"]).reset_index(drop=True)


def write_population_excel(out, df, ok, cs, covdf):
    key = ["Organism", "Chromosome", "Length_bp", "Feature", "Status", "N_plus", "N_minus", "Bins_used",
           "Len_balance_pct", "Len_KS_D", "Len_KS_D_expected", "Len_D_over_expected", "Len_r"]
    for k in PAIR_KEYS:
        key += [f"{k}_delta", f"{k}_r", f"{k}_delta_self", f"{k}_r_self", f"{k}_rank_d", f"{k}_rank_r",
                f"{k}_mean_plus", f"{k}_mean_minus"]
    for k in ("AT", "GC"):
        key += [f"{k}_dev_plus", f"{k}_dev_minus", f"{k}_dev_pooled", f"{k}_noise", f"{k}_silencing", f"{k}_r_opp"]
    key += ["Cross_r_nearest_of7", "Cross_d_nearest_of7", "Cross_r_nearest_of4", "Cross_d_nearest_of4"]
    key = [c for c in key if c in df.columns]
    readme = [f"Supplementary tables of the population analysis (Figure4_final v{VERSION})", "",
              "Sheet S2a summary by class: median, quartiles and range across chromosomes per class (the numbers for the Results text).",
              "Sheet S2b by organism: medians per organism and class.",
              "Sheet S2c chromosomes: one row per chromosome and class (Status = insufficient when < 50 blocks per strand or "
              "< 8 usable length bins).",
              "Sheet S2d coverage: FASTA records read, rejected and analysed per organism.", ""] + README_COMMON
    return write_excel(Path(out) / "Supp_Tables_population.xlsx",
                       {"S2a summary by class": cs, "S2b by organism": organism_summary(ok),
                        "S2c chromosomes": df[key], "S2d coverage": covdf}, readme)


def summarise_population(make_fig=True, extended=False):
    df, covdf = load_results()
    if df is None:
        print("[ERROR] no finished organisms in", org_dir())
        return
    out = Path(OUTPUT_DIR)
    out.mkdir(parents=True, exist_ok=True)
    ok = df[df["Status"] == "ok"].copy()
    ok["Feature"] = ok["Feature"].astype(str)
    if ok.empty:
        print("[ERROR] no usable chromosome x class rows")
        return
    df.to_csv(out / "per_chromosome_similarity.csv", index=False, float_format="%.6g")
    report, cs = build_population_report(df, ok, covdf)
    cs.to_csv(out / "summary_by_class.csv", index=False, float_format="%.6g")
    if covdf is not None and len(covdf):
        covdf.to_csv(out / "coverage_by_organism.csv", index=False)
    (out / "Population_report.txt").write_text(report, encoding="utf-8")
    print("\n" + report)
    written = ["per_chromosome_similarity.csv", "summary_by_class.csv", "Population_report.txt"]
    if make_fig:
        figs = [supp_s1_population_7col, supp_s2_population_matrix]
        if extended:
            figs += [supp_s1_length, supp_s2_pairs, supp_s3_cross, supp_s4_ratios]
        for fn in figs:
            written += [p.name for p in fn(ok, out)]
    if WRITE_EXCEL:
        x = write_population_excel(out, df, ok, cs, covdf)
        if x is not None:
            written.append(x.name)
    print(f"[DONE] written to {out}: " + ", ".join(written))


# =============================================================================
# 13. MAIN MODE (one chromosome -> Figure 4)
# =============================================================================

def run_main(no_figure, footnote):
    gdir = Path(INPUT_ROOT) / ORGANISM
    if not gdir.is_dir():
        print(f"[ERROR] organism folder not found: {gdir}")
        return 1
    out_dir = Path(OUTPUT_DIR) / ORGANISM
    out_dir.mkdir(parents=True, exist_ok=True)
    wanted = {MAIN_CHROMOSOME} if MAIN_CHROMOSOME else None
    params_line = (f"Length bins {BIN_WIDTH_NT} nt; >= {MIN_ELEMENTS_PER_STRAND_BIN} blocks per strand per bin; tail trimming "
                   f"{'on' if TRIM_LENGTH_RANGE else 'off'} (gap {MAX_GAP_BINS}); blocks > {MAX_FEATURE_LENGTH} nt dropped; "
                   f">= {MIN_USABLE_BINS} usable bins per class")
    print("=" * 90)
    print(f"[START] Figure 4 final v{VERSION}   organism={ORGANISM}   chromosome="
          f"{MAIN_CHROMOSOME or 'first valid chromosome'}")
    print(f"[INFO] {params_line}")
    print("[INFO] length windows (nt): " + "  ".join(f"{f} {class_window(f)[0]}-{class_window(f)[1]}" for f in FEATURES))
    print("=" * 90, flush=True)
    cov = defaultdict(int)
    for fpath in find_fasta_files(str(gdir)):
        for seq_id, header, blocks, length in iter_chromosomes(fpath, wanted, cov):
            print(f"[CHR] {seq_id}  {length:,} bp   ({header[:80]})", flush=True)
            res = {}
            for f in FEATURES:
                res[f] = analyse_feature(f, blocks[f][0], blocks[f][1], keep=True)
                print(f"  {f:<7} N+ {res[f]['n_p']:>7,}  N- {res[f]['n_m']:>7,}  "
                      f"{'' if res[f]['status'] == 'ok' else '(insufficient)'}", flush=True)
            blocks = None
            report, pairs, ratios, wide, mats = build_chr_report(res, seq_id, params_line)
            tag = safe_name(seq_id)
            wide.to_csv(out_dir / f"Figure4_stats_{tag}.csv", index=False, float_format="%.6g")
            pairs.to_csv(out_dir / f"Figure4_pairs_{tag}.csv", index=False, float_format="%.6g")
            ratios.to_csv(out_dir / f"Figure4_ratios_{tag}.csv", index=False, float_format="%.6g")
            mats.to_csv(out_dir / f"Figure4_cross_matrices_{tag}.csv", index=False, float_format="%.6g")
            (out_dir / f"Figure4_report_{tag}.txt").write_text(report, encoding="utf-8")
            print("\n" + report)
            written = [f"Figure4_stats_{tag}.csv", f"Figure4_pairs_{tag}.csv", f"Figure4_ratios_{tag}.csv",
                       f"Figure4_cross_matrices_{tag}.csv", f"Figure4_report_{tag}.txt"]
            if WRITE_EXCEL:
                x = write_chr_excel(out_dir, tag, seq_id, wide, pairs, ratios, cross_df_for_excel(res), mats)
                if x is not None:
                    written.append(x.name)
            if not no_figure:
                written += [p.name for p in draw_figure4(res, seq_id, out_dir, f"Figure_4_{ORGANISM}_{tag}", footnote)]
                written += [p.name for p in draw_pattern_matrices(res, seq_id, out_dir)]
            print(f"[DONE] written to {out_dir}: " + ", ".join(written))
            return 0
    print("[ERROR] chromosome not found (check --organism, --chromosome, --min-chromosome-bp)")
    return 1


def main():
    global INPUT_ROOT, OUTPUT_DIR, ORGANISM, MAIN_CHROMOSOME, MIN_CHROMOSOME_BP, CHUNK_MB, WORKERS
    global SAME_RANGE_FOR_ALL, TRIM_LENGTH_RANGE, LABELS_MODE, WRITE_EXCEL
    ap = argparse.ArgumentParser(description=__doc__.split("Usage")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=("main", "population", "summarise", "all"), default="all")
    ap.add_argument("--input-root", default=INPUT_ROOT)
    ap.add_argument("--output-dir", default=OUTPUT_DIR)
    ap.add_argument("--organism", default=ORGANISM, help="organism folder for --mode main")
    ap.add_argument("--chromosome", default=MAIN_CHROMOSOME,
                    help="sequence id for --mode main (default: first valid chromosome = chromosome 1)")
    ap.add_argument("--organisms", nargs="*", default=[], help="restrict --mode population to these folders")
    ap.add_argument("--min-chromosome-bp", type=int, default=MIN_CHROMOSOME_BP)
    ap.add_argument("--workers", type=int, default=WORKERS)
    ap.add_argument("--chunk-mb", type=float, default=CHUNK_MB)
    ap.add_argument("--fresh", action="store_true", help="redo organisms that are already finished")
    ap.add_argument("--no-figure", action="store_true")
    ap.add_argument("--footnote", default="")
    ap.add_argument("--same-range", nargs=2, type=int, metavar=("LO", "HI"), default=None,
                    help="ONE length range (nt) for every class, e.g. --same-range 20 1000")
    ap.add_argument("--no-trim", action="store_true",
                    help="do not cut the sparse tail: the length range is exactly the one set in the script / --same-range")
    ap.add_argument("--extended", action="store_true", help="population mode: also write the detailed figures S3-S6")
    ap.add_argument("--labels", choices=("short", "full"), default=LABELS_MODE,
                    help="Figure 4 panel titles: short = minimal statistics (default), full = all statistics")
    ap.add_argument("--no-excel", action="store_true", help="do not write the Supp_Tables_*.xlsx workbooks")
    a = ap.parse_args()
    LABELS_MODE = a.labels
    WRITE_EXCEL = not a.no_excel
    if a.same_range:
        SAME_RANGE_FOR_ALL = (a.same_range[0], a.same_range[1])
    if a.no_trim:
        TRIM_LENGTH_RANGE = False
    INPUT_ROOT, OUTPUT_DIR, ORGANISM, MAIN_CHROMOSOME = a.input_root, a.output_dir, a.organism, a.chromosome
    MIN_CHROMOSOME_BP, CHUNK_MB, WORKERS = a.min_chromosome_bp, a.chunk_mb, max(1, a.workers)
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    rc = 0
    if a.mode in ("main", "all"):
        rc = run_main(a.no_figure, a.footnote)
    if a.mode in ("population", "all"):
        print("=" * 90)
        print(f"[START] Population run v{VERSION}   input={INPUT_ROOT}   workers={WORKERS}")
        print("=" * 90, flush=True)
        run_population(WORKERS, a.fresh, a.organisms)
    if a.mode in ("population", "summarise", "all"):
        summarise_population(make_fig=not a.no_figure, extended=a.extended)
    return rc


if __name__ == "__main__":
    sys.exit(main())