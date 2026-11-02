# RepeatMasker — mammal pipeline (install and run)

This folder turns the downloaded mammal genomes into the two inputs that Figures 3, 4 and 5
need:

```
Fasta_GFF_inputs/Chordata/ncbi_dataset/data/<GCF_...>/   (from Download_Data)
        │  Automation.py  (copy + split + RepeatMasker + annotation)
        ▼
Fasta_GFF_inputs/Mammals/<Organism_Name>/                <chr>.fna  <chr>.fna.out  <chr>.csv  *.gff.gz
        │                                   │
        │ transfer_CSV_files.py             │ ../Masking_GFF_OUT/produce_masked_fna.py
        ▼                                   ▼
Dataframes/Mammals_CSV/<organism>/<chr>.csv        Fasta_GFF_inputs/Masked_Fasta/<organism>/<chr>.fna.gz
        (Figure 3)                                        (Figures 4, 5, 8)
```

## 1. Install RepeatMasker (Ubuntu)

Follow the official instructions (https://github.com/Dfam-consortium/RepeatMasker) — this
is a summary of what they require:

| Needed | Notes |
|---|---|
| Perl ≥ 5.8.0, Python ≥ 3 | `sudo apt install perl python3 python3-h5py` (h5py is used to read the Dfam library files) |
| TRF ≥ 4.09 (Tandem Repeats Finder) | download the Linux binary from its GitHub releases (https://github.com/Benson-Genomics-Lab/TRF), rename it `trf`, `chmod +x`, put it on the `PATH` |
| a search engine | **RMBlast** (`rmblastn`), which this pipeline uses (`-engine rmblast`); download from https://www.repeatmasker.org/rmblast/ |
| repeat library | RepeatMasker has none built in: install **FamDB / Dfam partitions** (https://github.com/Dfam-consortium/FamDB/releases) into `RepeatMasker/Libraries/famdb/` — the manuscript used Dfam 3.8 — or give your own library with `-lib` |

Then, in the RepeatMasker folder (do not unpack into a folder that already contains a
`RepeatMasker` directory):

```bash
cd /usr/local                                  # or ~/local_tools
tar xzf RepeatMasker-<version>.tar.gz
cd RepeatMasker
perl ./configure                               # asks for the paths of trf and rmblastn
./RepeatMasker -help                           # test
```

Tell this repository where RepeatMasker is (only needed if it is not on the `PATH`):

```bash
export CSPR_REPEATMASKER=/usr/local/RepeatMasker/RepeatMasker
export CSPR_RM_THREADS=8          # value of RepeatMasker -pa  (default 8)
```

> `-pa N` is RepeatMasker's parallel setting; choose a value that fits your machine
> (see `RepeatMasker -help`).

## 2. Everything in one go

Prerequisite: the mammals are downloaded (`Download_Data/README.md`) into
`Fasta_GFF_inputs/Chordata/`.

```bash
python3 RepeatMasker/Automation.py                      # all organisms of the list
python3 RepeatMasker/Automation.py --only Mus_musculus  # one organism, good for a first test
python3 RepeatMasker/Automation.py --skip-annotation    # only up to RepeatMasker
```

The organism list is `RepeatMasker/mamals_refrence_annotated_2025_7_12.csv` (columns
*Organism Name*, *Organism Taxonomic ID*, assembly accession in the first column); use
`--list` for another file (`.csv` or `.tsv`). The run can be stopped and started again.

Step 4 of `Automation.py` (per-chromosome annotation → `<chr>.csv`) needs
`RepeatMasker/eukaryotes_for_chargaff.py`. If that file is missing the script stops with a
message at that point; use `--skip-annotation` until it is added.

## 3. After the run

```bash
python3 RepeatMasker/transfer_CSV_files.py              # CSV tables -> Dataframes/Mammals_CSV  (Figure 3)
python3 Masking_GFF_OUT/produce_masked_fna.py           # masked FASTA  -> Fasta_GFF_inputs/Masked_Fasta  (Figures 4, 5)
python3 RepeatMasker/mix_repeatmasker_result.py         # optional: which organisms are unfinished
python3 RepeatMasker/remove_some_results.py             # optional: dry run of the clean-up (add --delete to delete)
```

## 4. The single steps

| Script | What it does |
|---|---|
| `Genome_file_dividing.py <organism>` | splits `<genome>.fna(.gz)` of `Fasta_GFF_inputs/Mammals/<organism>/` into one `<chr>.fna` per chromosome in the same folder (scaffolds and patches are skipped; `--delete-original` removes the whole-genome file) |
| `Linux_in_python_RepeatMasker.py <organism> <taxon>` | runs `RepeatMasker -species <taxon> -no_is -nolow -norna -pa N -engine rmblast` on each `<chr>.fna` → `<chr>.fna.out`; chromosomes that already have a `.out` are skipped |
| `linux_in_python_with_espicial_lib.py <organism> --lib <file>` | same with a custom library (`-lib`) instead of `-species` |
| `filter_repeatmasker_library.py <fasta>` | lists the taxon names in a Dfam FASTA library (for building a custom library) |
| `mix_repeatmasker_result.py` | table of organisms with a different number of `.fna` and `.out` files |
| `transfer_CSV_files.py` | copies the `<chr>.csv` tables to `Dataframes/Mammals_CSV/` |
| `remove_some_results.py` | **dry run by default**; deletes `*.cat.gz`, `*.masked` and `*.csv` only with `--delete` (run `transfer_CSV_files.py` first!) |

`-species` is given the NCBI taxon ID (as in the original runs). If your RepeatMasker
version does not accept an ID, use the species name instead.
