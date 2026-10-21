# Dataframes/ — input tables (.csv)

Put the data tables here. They are **not** stored in git (see `.gitignore`);
download them from the data archive (Zenodo) or copy them from your own analysis.
The scripts find them through `config.py`, so nothing else needs editing.

| File / folder | Used by |
|---|---|
| `Cleaned_datafram_whole_organims_2026_07.csv` (master table, one row per chromosome) | Figures 1, 2, 6, 7 |
| `Viruses_2026-01-26_16-51-15.csv` | Figure 2 |
| `Refrence_viruses_list.csv` | Figure 2 |
| `Plasmides_2026-02-05_23-22-23_add_columns.csv` | Figure 2 |
| `simulation_results.csv` (artificial-chromosome simulation) | Figure 2 |
| `Mammals/<organism>/<chromosome>.csv` (per-chromosome element composition) | Figure 3 |
| `figure8_equilibrium_other_groups/<Group>_sampling_results.csv` (Bacteria, Archaea, Fungi, Arthropoda, Plantae, Chordata, Mammals) | Figure 8 |
| `fungi_ref_annotated_2025_10_14.tsv` | `GFF_Mapping/Fungi/new_fungal_CDS_UTR_Intron.py` |

The folder can live somewhere else: set the environment variable `CSPR_DATAFRAMES`
(or `CSPR_MASTER_CSV` for just the master table). See the top-level README.
