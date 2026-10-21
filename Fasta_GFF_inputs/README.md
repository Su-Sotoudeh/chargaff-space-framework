# Fasta_GFF_inputs/ — genome files (FASTA / GFF)

Put the genome files here. They are **not** stored in git (see `.gitignore`).

```
Fasta_GFF_inputs/
├── Masked_Fastas/                      RepeatMasker-masked, symbol-encoded FASTA
│   └── <organism>/                     one folder per organism
│       └── <chromosome>.fna(.gz) ...   (.fna / .fa / .fasta, optionally gzipped)
│
└── <Group>/ncbi_dataset/data/          unzipped NCBI genome package, one folder per group
    └── <assembly accession>/
        ├── *.fna                       genome sequence
        └── *.gff                       annotation
```

`<Group>` is one of: `Bacteria`, `Archaea`, `Fungi`, `Chordata`, `Arthropoda`,
`Plantae`, `Mitochondrion`, `Chloroplast`.
For `Mitochondrion` and `Chloroplast` the file `organelle_gene_report.csv` also sits in
`ncbi_dataset/data/`.

An NCBI package downloaded with `datasets download genome ...` and unzipped into
`Fasta_GFF_inputs/<Group>/` already has the `ncbi_dataset/data/...` layout.

| Folder | Used by |
|---|---|
| `Masked_Fastas/` | Figure 4, Figure 5 (`count_from_masked_file.py`, `Figure_5_<class>.py`), Figure 8 (`Genes_infor_extraction_Mammals.py`), `GFF_Mapping/Mammals/` |
| `<Group>/ncbi_dataset/data/` | `GFF_Mapping/<Group>/`, Figure 8 (`Genes_info_extraction.py`, set `NAME` at the top) |

The folder can live somewhere else (e.g. a large drive): set the environment variable
`CSPR_FASTA_GFF`. See the top-level README.
