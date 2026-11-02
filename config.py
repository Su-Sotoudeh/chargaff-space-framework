"""
config.py -- the ONE place that knows where files live.

Every script in this repository imports its input and output locations from
here; there are no personal / absolute paths anywhere else.  Everything is
computed relative to this file, so the repository works wherever it is cloned
(laptop, lab PC, HPC) on Linux, macOS and Windows.

Folder layout (relative to the repository root)
-----------------------------------------------
    Dataframes/                 INPUT  tables (.csv) used by the figure scripts
        Cleaned_datafram_whole_organims_2026_07.csv   <- master table (Figures 1,2,6,7)
        simulation_results.csv                        <- Figure 2
        Mammals_CSV/<organism>/<chromosome>.csv       <- Figure 3
        figure8_equilibrium_other_groups/*_sampling_results.csv  <- Figure 8

    Fasta_GFF_inputs/           INPUT  genome files (FASTA / GFF)
        Mammals/<organism>/<chromosome>.fna(.gz)      <- RepeatMasker work folder (whole genome -> per-chromosome files)
        Masked_Fasta/<organism>/<chromosome>.fna.gz   <- Figures 4, 5, 8, GFF_Mapping/Mammals
        <Group>/ncbi_dataset/data/<assembly>/*.fna|*.gff   <- GFF_Mapping/<Group>, Figure 8
            where <Group> is Bacteria, Archaea, Fungi, Chordata, Arthropoda,
            Plantae, Mitochondrion, Chloroplast, ...  (an unzipped NCBI
            `datasets download genome ...` package keeps this exact layout)

    Figure_<n>/outputs/         OUTPUT everything a figure script writes
    GFF_Mapping/outputs/<Group>/ OUTPUT tables written by the GFF mapping scripts
    (all `outputs/` folders are git-ignored and are re-created automatically)

Overriding locations (optional)
-------------------------------
Large data does not have to live inside the repository.  Set any of these
environment variables to point somewhere else, e.g. on the HPC:

    CSPR_DATAFRAMES   folder with the input tables      (default: <repo>/Dataframes)
    CSPR_FASTA_GFF    folder with the genome files      (default: <repo>/Fasta_GFF_inputs)
    CSPR_OUTPUTS      one folder to receive ALL outputs (default: <repo>/<Figure_n>/outputs)
    CSPR_MASTER_CSV   full path of the master table     (default: <Dataframes>/Cleaned_datafram_whole_organims_2026_07.csv)
    CSPR_REPEATMASKER RepeatMasker program (default: `RepeatMasker` on PATH)
    CSPR_RM_THREADS   RepeatMasker -pa value            (default: 8)
    CHARGAFF_CSV      same as CSPR_MASTER_CSV (kept for older command lines)

    Linux/macOS:  export CSPR_FASTA_GFF=/data/genomes
    Windows:      set CSPR_FASTA_GFF=D:\\genomes
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def _env_path(name: str, default: Path) -> Path:
    value = os.environ.get(name)
    return Path(value).expanduser() if value else Path(default)


# ----------------------------------------------------------------------------
# INPUT locations
# ----------------------------------------------------------------------------
DATAFRAMES = _env_path("CSPR_DATAFRAMES", ROOT / "Dataframes")
FASTA_GFF = _env_path("CSPR_FASTA_GFF", ROOT / "Fasta_GFF_inputs")

# Master table (one row per chromosome) -- Figures 1, 2, 6, 7
_master_default = DATAFRAMES / "Cleaned_datafram_whole_organims_2026_07.csv"
MASTER_CSV = Path(
    os.environ.get("CSPR_MASTER_CSV") or os.environ.get("CHARGAFF_CSV") or _master_default
).expanduser()

# Simulation table -- Figure 2
SIMULATION_CSV = DATAFRAMES / "simulation_results.csv"

# Per-chromosome composition tables (Mammals_CSV/<organism>/<chromosome>.csv) -- Figure 3
MAMMALS_CHR_DIR = DATAFRAMES / "Mammals_CSV"

# Sampling-experiment tables (<Group>_sampling_results.csv) -- Figure 8
FIG8_SAMPLING_DIR = DATAFRAMES / "figure8_equilibrium_other_groups"

# NCBI reference lists used by the Fungi mapping script
FUNGI_REFERENCE_LIST = DATAFRAMES / "fungi_ref_annotated_2025_10_14.tsv"

# RepeatMasker-masked, symbol-encoded FASTA files (one sub-folder per organism)
MASKED_FASTAS = FASTA_GFF / "Masked_Fasta"

# Mammal work folder: whole-genome FASTA is split into one file per chromosome
# here (RepeatMasker/Genome_file_dividing.py), RepeatMasker is run on these files,
# and the .fna.out results stay next to them.
MAMMALS_WORK_DIR = FASTA_GFF / "Mammals"

# RepeatMasker installation (see RepeatMasker/README.md).  Override with env vars.
REPEATMASKER_BIN = os.environ.get("CSPR_REPEATMASKER", "RepeatMasker")   # on PATH, or a full path
REPEATMASKER_THREADS = int(os.environ.get("CSPR_RM_THREADS", "8"))         # RepeatMasker -pa value


def ncbi_data_dir(group: str) -> Path:
    """Folder that holds the per-assembly sub-folders of an NCBI genome package.

    Fasta_GFF_inputs/<group>/ncbi_dataset/data/<assembly>/{*.fna, *.gff}
    """
    return FASTA_GFF / group / "ncbi_dataset" / "data"


# ----------------------------------------------------------------------------
# OUTPUT locations
# ----------------------------------------------------------------------------
def figure_dir(n: int) -> Path:
    """Folder of Figure_n."""
    return ROOT / f"Figure_{n}"


def output_dir(key, *subfolders: str) -> Path:
    """Output folder (created automatically).

    output_dir(4)                       -> Figure_4/outputs
    output_dir(8, "ML")                 -> Figure_8/outputs/ML
    output_dir("GFF_Mapping", "Fungi")  -> GFF_Mapping/outputs/Fungi

    If CSPR_OUTPUTS is set, everything goes below that folder instead:
    $CSPR_OUTPUTS/Figure_4, $CSPR_OUTPUTS/Figure_8/ML, ...
    """
    name = f"Figure_{key}" if isinstance(key, int) else str(key)
    override = os.environ.get("CSPR_OUTPUTS")
    base = (Path(override).expanduser() / name) if override else (ROOT / name / "outputs")
    d = base.joinpath(*subfolders)
    d.mkdir(parents=True, exist_ok=True)
    return d


def figure_input(n: int, filename: str) -> Path:
    """A table a figure reads that may be produced by an earlier step OR shipped with the repository.

    Looked up first in the figure's output folder (what the pipeline wrote), then in the figure's
    own folder (the small test file committed to GitHub).  Stops with a clear message if neither exists.
    """
    for d in (output_dir(n), figure_dir(n)):
        if (d / filename).exists():
            return d / filename
    raise SystemExit(
        f"\n[config] Input not found: {filename}\n"
        f"         Looked in {output_dir(n)} and {figure_dir(n)}.\n"
        f"         Run the Figure_{n} counting scripts first, or put the file in {figure_dir(n)}.\n")


def require(path, hint: str = "") -> Path:
    """Return `path` if it exists, otherwise stop with a message saying what to put where."""
    p = Path(path)
    if not p.exists():
        raise SystemExit(
            f"\n[config] Input not found: {p}\n"
            f"         {hint or 'See README.md ("Where to put your data").'}\n"
        )
    return p
