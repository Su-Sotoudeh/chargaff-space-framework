# Data download guide

This guide explains how to download the genome data used in the Chargaff Space Framework analyses: **cellular genomes** (six taxon groups), and **mitochondria, chloroplasts (plastids), plasmids and viruses**. Everything is saved under `Fasta_GFF_inputs/` (see `paths.py`).

> **Why dates matter.** NCBI updates its databases continuously, so the number of genomes changes from month to month. Every download is stamped with the date and the `datasets` version in `download_log.tsv` and in `MANIFEST.txt`, and the data used for the paper should be archived (Zenodo) rather than re-downloaded.

## 1. Overview

| Data | Source | Method | Script |
|---|---|---|---|
| Cellular genomes (Chordata, Viridiplantae, Arthropoda, Fungi, Bacteria, Archaea) | NCBI Datasets (RefSeq reference genomes, annotated, complete or chromosome level) | dehydrate → unzip → rehydrate, in a loop | `scripts/download/download_cellular.sh` |
| Mitochondria | NCBI RefSeq release (`mitochondrion`) | FTP, reads the server listing | `scripts/download/download_organelles_plasmids_viruses.sh` |
| Chloroplasts / plastids | NCBI RefSeq release (`plastid`) | same | same script |
| Plasmids | NCBI RefSeq `genomes/refseq/plasmid` | same | same script |
| Viruses | NCBI RefSeq release (`viral`) | same | same script |
| Virus genome type and reference taxonomy | ICTV Master Species List | manual download (spreadsheet) | none |

## 2. Requirements

- Linux or macOS shell (on Windows use WSL or Git Bash).
- [NCBI Datasets command-line tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) (`datasets`), plus `unzip`, `curl` and `gzip`.
- Disk space: the cellular genomes need hundreds of GB. Check free space first with `df -h .`.
- A stable connection. All scripts can be re-run safely and resume interrupted work.

Check the installation:

```bash
datasets --version
unzip -v | head -n1
curl --version | head -n1
```

## 3. Cellular genomes (loop over taxon groups)

### What each step does

Large genome sets are downloaded in three steps. The *dehydrated* package contains only the metadata and the list of files, so it is small and fast; the real sequence files are fetched afterwards in parallel, and an interrupted download can be resumed.

| Step | Purpose | Command (one line) |
|---|---|---|
| 1. Dehydrated download | get the package with the file list | `datasets download genome taxon "Archaea" --reference --annotated --assembly-level complete,chromosome --include genome,gff3 --dehydrated --filename Archaea_2026_11.zip` |
| 2. Unzip | extract the package | `unzip Archaea_2026_11.zip -d Archaea_2026_11` |
| 3. Rehydrate | download FASTA and GFF3 files | `datasets rehydrate --directory Archaea_2026_11 --max-workers 8` |

Options used:

- `--reference` keeps one reference genome per species.
- `--annotated` keeps only annotated assemblies (needed for the +annotation / −annotation analysis).
- `--assembly-level complete,chromosome` excludes scaffold- and contig-level assemblies.
- `--include genome,gff3` downloads the sequence and the annotation.

> Write each command on **one line** (or use a trailing backslash `\` with nothing after it). A line break in the middle of a command is the most common cause of errors when copying.

### Run all groups with the loop script

```bash
bash scripts/download/download_cellular.sh
```

This loops over `Chordata Viridiplantae Arthropoda Fungi Bacteria Archaea` and runs the three steps for each one. To download only some groups, or a different group, list them:

```bash
bash scripts/download/download_cellular.sh Archaea Fungi
```

Useful settings (environment variables):

```bash
STAMP=2026_11 WORKERS=16 bash scripts/download/download_cellular.sh
```

What the script does for you:

- Creates `Fasta_GFF_inputs/cellular/<Group>_<YYYY_MM>/`.
- Skips groups that are already complete, and resumes interrupted ones.
- Continues with the next group if one fails, and lists the failed groups at the end.
- Appends a line per group to `Fasta_GFF_inputs/cellular/download_log.tsv` (date, group, `datasets` version, number of assemblies).

If you want the same loop inline, without the script:

```bash
for TAXON in Chordata Viridiplantae Arthropoda Fungi Bacteria Archaea; do
  NAME="${TAXON}_$(date +%Y_%m)"
  datasets download genome taxon "$TAXON" --reference --annotated --assembly-level complete,chromosome --include genome,gff3 --dehydrated --filename "$NAME.zip"
  unzip -q "$NAME.zip" -d "$NAME"
  datasets rehydrate --directory "$NAME" --max-workers 8
done
```

Mammals are part of Chordata. Separate the mammalian subset afterwards from the taxonomy information in the downloaded metadata.

## 4. Mitochondria and chloroplasts (plastids)

NCBI's organelle web page (`Datasets → Organelle`) is designed for browsing and manual download. For a reproducible command-line download, the script uses the RefSeq release folders on the NCBI FTP server instead:

```bash
bash scripts/download/download_organelles_plasmids_viruses.sh mitochondrion plastid
```

Results go to `Fasta_GFF_inputs/mitochondrion/` and `Fasta_GFF_inputs/plastid/`:

- `*.genomic.fna.gz`: sequences
- `*.genomic.gbff.gz`: GenBank flat files with the annotation (genes, strand)
- `MANIFEST.txt`: file list and SHA-256 checksums with the download date

Notes:

- These folders contain **all** RefSeq organelle genomes (animals, fungi, plants, protists). Filter to the taxa you need (for example eukaryotes, taxon 2759) in the analysis from the organism information in the GenBank records.
- `plastid` includes chloroplasts and other plastids. Select chloroplasts in the analysis from the record description.
- The release files provide annotation in GenBank format. If your scripts need GFF3, convert the GenBank files, or download specific assemblies as GFF3 from the NCBI web page.
- Web alternatives, if you prefer to download by hand:
  - Mitochondria: https://www.ncbi.nlm.nih.gov/datasets/organelle/?taxon=2759&type=Mitochondrion
  - Chloroplasts: https://www.ncbi.nlm.nih.gov/datasets/organelle/?taxon=2759&type=Chloroplast

## 5. Plasmids

```bash
bash scripts/download/download_organelles_plasmids_viruses.sh plasmid
```

Source: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/. Output: `Fasta_GFF_inputs/plasmid/` (FASTA and GenBank files plus `MANIFEST.txt`).

## 6. Viruses

```bash
bash scripts/download/download_organelles_plasmids_viruses.sh viral
```

Source: https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/. Output: `Fasta_GFF_inputs/viral/`.

The earlier web query restricted the list to annotated RefSeq/GenBank viral assemblies at complete or chromosome level (https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=10239&annotated_only=true&refseq_annotation=true&genbank_annotation=true&assembly_level=2:3). The RefSeq release folder is a reproducible equivalent for RefSeq records; apply the annotation and completeness filters in the analysis.

### Virus genome type and reference taxonomy (ICTV)

The genome type (for example dsDNA, ssRNA(+)) comes from the ICTV Master Species List:

1. Open https://ictv.global/msl.
2. Download the latest Master Species List spreadsheet.
3. Save it as `Fasta_GFF_inputs/viral/ICTV/ICTV_MSL_<release>.xlsx` and note the release number in the README of that folder.

## 7. Run all four organelle, plasmid and viral sets at once

```bash
bash scripts/download/download_organelles_plasmids_viruses.sh
```

## 8. Check the downloads

```bash
bash scripts/download/check_downloads.sh
```

This reports, for each cellular group, the number of assemblies and how many have FASTA and GFF files (assemblies without GFF are listed in `<Group>_missing_gff.txt`), and tests every `.gz` file of the other sets for corruption.

## 9. Expected folder structure

```
Fasta_GFF_inputs/
├── cellular/
│   ├── download_log.tsv
│   ├── Archaea_2026_11/ncbi_dataset/data/GCF_.../{*.fna, genomic.gff}
│   └── ... (one folder per group)
├── mitochondrion/   *.genomic.fna.gz, *.genomic.gbff.gz, MANIFEST.txt
├── plastid/         same
├── plasmid/         same
└── viral/           same, plus ICTV/ICTV_MSL_<release>.xlsx
```

## 10. Troubleshooting

| Problem | What to do |
|---|---|
| `datasets: command not found` | Install the NCBI Datasets CLI and make sure it is on your `PATH`. |
| Download stops or the connection drops | Run the same command again; finished files are kept and the rest is resumed. |
| `--max-workers` not accepted | Run `datasets rehydrate --help` and use the option name shown there, or remove the option. |
| `no ... files found; the server layout may have changed` | Open the URL in a browser, check the new file names, and adjust the pattern in the script. |
| A group fails repeatedly | Delete the `<Group>_<date>` folder and its `.zip`, then run the script for that group again. |
| Disk full | Free space or set `CSF_FASTA_GFF` to a larger disk. |

## 11. Reproducibility and archiving

- Keep `download_log.tsv` and each `MANIFEST.txt`; they record dates, software version and checksums.
- State the download month in the manuscript's data availability section.
- Archive the exact data (or the processed tables derived from it) on Zenodo or Figshare, and link it from the main README.
