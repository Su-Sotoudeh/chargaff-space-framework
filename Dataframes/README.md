# Dataframes/ — input tables (.csv)

Input tables of the figure scripts. The master table, `simulation_results.csv` and
`Mammals_CSV/` are small and are committed to git; the heavy Figure 8 tables in
`figure8_equilibrium_other_groups/` are not (produce them with the Figure 8 scripts, or
download them from the data archive).
The scripts find them through `config.py`, so nothing else needs editing.

| File / folder | Used by |
|---|---|
| `Cleaned_datafram_whole_organims_2026_07.csv` (master table, one row per chromosome; only the columns the scripts use) | Figures 1, 2, 6, 7 |
| `simulation_results.csv` (artificial-chromosome simulation) | Figure 2 |
| `Mammals_CSV/<organism>/<chromosome>.csv` (per-chromosome element composition; made by `RepeatMasker/transfer_CSV_files.py`) | Figure 3 |
| `figure8_equilibrium_other_groups/<Group>_sampling_results.csv` (Bacteria, Archaea, Fungi, Arthropoda, Plantae, Chordata, Mammals) | Figure 8 |
| `fungi_ref_annotated_2025_10_14.tsv` | `GFF_Mapping/Fungi/new_fungal_CDS_UTR_Intron.py` |

The folder can live somewhere else: set the environment variable `CSPR_DATAFRAMES`
(or `CSPR_MASTER_CSV` for just the master table). See the top-level README.
