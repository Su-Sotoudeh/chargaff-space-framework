# Fasta_GFF_inputs/ — genome files (FASTA / GFF)

Put the genome files here. They are **not** stored in git (see `.gitignore`).

```
Fasta_GFF_inputs/
├── Mammals/<organism>/                 RepeatMasker work folder: <chr>.fna, <chr>.fna.out, <chr>.csv, *.gff.gz
│                                       (made by RepeatMasker/Automation.py)
│
├── Masked_Fasta/                       RepeatMasker-masked, symbol-encoded FASTA
│   └── <organism>/<chromosome>.fna.gz  one folder per organism, one file per chromosome
│                                       (made by Masking_GFF_OUT/produce_masked_fna.py)
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

Download everything with the scripts in `Download_Data/` (see `Download_Data/README.md`);
they create exactly this layout. Mammal genomes are the part of `Chordata` listed in
`RepeatMasker/mamals_refrence_annotated_2025_7_12.csv`; see `RepeatMasker/README.md`.

| Folder | Used by |
|---|---|
| `Masked_Fasta/` | Figure 4, Figure 5 (`count_from_masked_file.py`, `Figure_5_<class>.py`), Figure 8 (`Genes_infor_extraction_Mammals.py`), `GFF_Mapping/Mammals/` |
| `Mammals/` | `RepeatMasker/`, `Masking_GFF_OUT/produce_masked_fna.py` |
| `<Group>/ncbi_dataset/data/` | `GFF_Mapping/<Group>/`, Figure 8 (`Genes_info_extraction.py`, set `NAME` at the top) |

The folder can live somewhere else (e.g. a large drive): set the environment variable
`CSPR_FASTA_GFF`. See the top-level README.
