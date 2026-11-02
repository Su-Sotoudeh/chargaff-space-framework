# Download_Data — getting the genomes from NCBI

Everything is saved under `Fasta_GFF_inputs/` (or the folder in `CSPR_FASTA_GFF`, see
`config.py`), in exactly the layout the analysis scripts read:

```
Fasta_GFF_inputs/
├── Chordata/ncbi_dataset/data/<GCF_...>/{*.fna, genomic.gff}      (mammals are part of Chordata)
├── Plantae/    Arthropoda/    Fungi/    Bacteria/    Archaea/        same layout
├── Mitochondrion/ncbi_dataset/data/...                             manual web download, section 4
└── Chloroplast/ncbi_dataset/data/...                               manual web download, section 4
```

> NCBI changes every month. Keep `Fasta_GFF_inputs/download_log.tsv` (date, `datasets`
> version, number of assemblies per group) and state the download month in the paper.

## 1. Install NCBI Datasets on Ubuntu

The tool is two programs: `datasets` (download) and `dataformat` (convert metadata).

**Option A — plain download (no conda):**

```bash
mkdir -p ~/.local/bin && cd ~/.local/bin
curl -o datasets   'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat
```

Ubuntu puts `~/.local/bin` on the `PATH` automatically after you log in again; for the
current terminal run `export PATH="$HOME/.local/bin:$PATH"`.

**Option B — conda:**

```bash
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli
```

Also needed: `sudo apt install unzip curl`.

Check:

```bash
datasets --version
dataformat --version
```

Official page: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/

## 2. Download the cellular genomes (one command)

```bash
bash Download_Data/scripts/download/download_cellular.sh                  # all six groups
bash Download_Data/scripts/download/download_cellular.sh Archaea Fungi   # only some
WORKERS=16 bash Download_Data/scripts/download/download_cellular.sh      # more parallel downloads
```

Run it from the repository root (or anywhere; it finds the repository itself). Groups:
`Chordata Viridiplantae Arthropoda Fungi Bacteria Archaea` — the NCBI taxon *Viridiplantae*
is stored in the folder **`Plantae`**.

What the script does for each group (the NCBI "dehydrate → rehydrate" workflow, which
makes big downloads resumable):

| Step | Command it runs |
|---|---|
| 1. dehydrated package (metadata and file list only, small) | `datasets download genome taxon <Group> --reference --annotated --assembly-level complete,chromosome --include genome,gff3 --dehydrated --filename <zip>` |
| 2. unzip into `Fasta_GFF_inputs/<Group>/` | `unzip` |
| 3. fetch the real FASTA and GFF3 files | `datasets rehydrate --directory Fasta_GFF_inputs/<Group> --max-workers 8` |

Options: `--reference` one reference genome per species; `--annotated` only annotated
assemblies (needed for the +/− annotation analysis); `--assembly-level complete,chromosome`
excludes scaffold/contig assemblies; `--include genome,gff3` sequence and annotation.

It is safe to re-run: finished groups are skipped, interrupted ones resume, a failed group
does not stop the others. It needs hundreds of GB of disk — check with `df -h .` and, if
needed, point `CSPR_FASTA_GFF` to a bigger drive.

## 3. Check the download

```bash
bash Download_Data/scripts/download/check_downloads.sh
```

Lists, per group, how many assemblies have a FASTA and a GFF file (assemblies without GFF
are written to `Fasta_GFF_inputs/<Group>_missing_gff.txt`).

## 4. Mitochondria and chloroplasts (manual, from the NCBI web page)

The command-line `datasets download` has no organelle command (only gene, genome, taxonomy,
virus), so organelles come from the web page:

1. Mitochondria: https://www.ncbi.nlm.nih.gov/datasets/organelle/?taxon=2759&type=Mitochondrion
   Chloroplasts: https://www.ncbi.nlm.nih.gov/datasets/organelle/?taxon=2759&type=Chloroplast
2. Use *Download* → choose genomic sequence (FASTA), annotation (GFF3) and the gene report.
3. Unzip into `Fasta_GFF_inputs/Mitochondrion/` (respectively `Chloroplast/`) so that
   `Fasta_GFF_inputs/Mitochondrion/ncbi_dataset/data/` exists. The scripts in
   `GFF_Mapping/Mito/` and `GFF_Mapping/Chrloplast/` also read
   `ncbi_dataset/data/organelle_gene_report.csv` from that folder.

## 5. Optional: RefSeq FTP files (plasmids, viruses, organelles in GenBank format)

```bash
bash Download_Data/scripts/download/download_organelles_plasmids_viruses.sh            # all four sets
bash Download_Data/scripts/download/download_organelles_plasmids_viruses.sh viral      # one set
```

Output: `Fasta_GFF_inputs/RefSeq_FTP/<set>/*.genomic.{fna,gbff}.gz` plus `MANIFEST.txt`
(checksums and date). These are *not* read by `GFF_Mapping/` (they are GenBank format).
Virus genome types come from the ICTV Master Species List (https://ictv.global/msl).

## 6. Troubleshooting

| Problem | What to do |
|---|---|
| `datasets: command not found` | `export PATH="$HOME/.local/bin:$PATH"` or activate the conda environment. |
| Download stops or the connection drops | Run the same command again; finished files are kept. |
| `--max-workers` not accepted | `datasets rehydrate --help` and use the option name shown there. |
| A group fails repeatedly | Delete `Fasta_GFF_inputs/<Group>/` and `Fasta_GFF_inputs/_dehydrated_zips/<Group>.zip`, run again. |
| Disk full | Free space, or set `CSPR_FASTA_GFF` to a larger disk. |
