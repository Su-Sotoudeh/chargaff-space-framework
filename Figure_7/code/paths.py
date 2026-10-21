"""
paths.py
---------
Output / input locations for the Figure 7 pipeline.  Everything is resolved
through the repository-level ``config.py`` so the pipeline works wherever the
repository is cloned (laptop, lab PC, HPC).

    input  : config.MASTER_CSV      (Dataframes/Cleaned_datafram_whole_organims_2026_07.csv)
    tables : Figure_7/outputs/tables/   (CSVs: scores, cleaned data, summaries)
    figures: Figure_7/outputs/figures/  (PNG + PDF figures)

Set CSPR_OUTPUTS to send all outputs to another drive (see config.py).
"""

import os
import sys
from pathlib import Path

CODE_DIR = os.path.dirname(os.path.abspath(__file__))

# make the repo-level config.py importable (code/ -> Figure_7/ -> repo root)
sys.path.insert(0, str(next(p for p in Path(CODE_DIR).parents if (p / "config.py").exists())))
import config  # noqa: E402

MASTER_CSV = str(config.MASTER_CSV)
TABLES_DIR = str(config.output_dir(7, "tables"))
FIGURES_DIR = str(config.output_dir(7, "figures"))
PROJECT_ROOT = str(config.figure_dir(7))


def table_path(filename):
    return os.path.join(TABLES_DIR, filename)


def figure_path(filename_no_ext):
    return os.path.join(FIGURES_DIR, filename_no_ext)
