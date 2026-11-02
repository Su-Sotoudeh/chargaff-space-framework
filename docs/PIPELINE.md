# Pipeline overview — from NCBI to the figures

```
NCBI ──Download_Data──► Fasta_GFF_inputs/<Group>/ncbi_dataset/data/…        (FASTA + GFF)
                              │
        ┌─────────────────────┼───────────────────────────────┐
        │ GFF_Mapping/<Group>/ (per-group tables)             │ mammals only
        ▼                                                     ▼
 master table (Dataframes/Cleaned_datafram_whole_organims_2026_07.csv)   RepeatMasker/Automation.py
        │                                                     │  → Fasta_GFF_inputs/Mammals/<org>/  <chr>.fna .out .csv
        │                                          ┌──────────┴───────────┐
        │                              transfer_CSV_files.py     Masking_GFF_OUT/produce_masked_fna.py
        │                                          ▼                       ▼
        │                              Dataframes/Mammals_CSV      Fasta_GFF_inputs/Masked_Fasta
        ▼                                          ▼                       ▼
 Figures 1, 2, 6, 7                          Figure 3              Figures 4, 5, 8
```

| Figure | Input | Run |
|---|---|---|
| 1, 6 | master table | `Figure_1/Figure_1.py`, `Figure_6/Figure_6.py` |
| 2 | master table + `simulation_results.csv` | `Figure_2/Figure_2_Final.py` |
| 3 | `Dataframes/Mammals_CSV/` | `Figure_3/Figure_3.py` |
| 4 | `Fasta_GFF_inputs/Masked_Fasta/<organism>/` | `Figure_4/Figure_4.py --organism <name>` |
| 5 A–C | `Figure_5/4k_kmers_normalised_MCRRC_PC0UW.csv` (committed test file) | `Figure_5/Panel_A.py`, `Panel_B.py`, `Panel_C.py` |
| 5 D | the same file + per-class tables | run `Figure_5/Figure_5_<CDS,UTR,INT,SINE,LINE,LTR,DNA>.py` (from `Masked_Fasta`), then `Panel_D.py` |
| 7 | master table | `Figure_7/code/run_pipeline.py` |
| 8 | `Dataframes/figure8_equilibrium_other_groups/*_sampling_results.csv` | `Figure_8/Figure_8_Panel_ABC_std.py` (the heavy tables are made by the other Figure 8 scripts) |

All scripts read and write through `config.py`; outputs go to `Figure_<n>/outputs/`
(Figure 5 panels and the per-class tables included).
