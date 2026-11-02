#!/usr/bin/env bash
# Download reference, annotated genomes (FASTA + GFF3) for the cellular taxon groups
# from NCBI with the dehydrate -> unzip -> rehydrate workflow.
#
# Result (the layout GFF_Mapping/ and config.py expect):
#   Fasta_GFF_inputs/<Group>/ncbi_dataset/data/<assembly accession>/{*.fna, genomic.gff}
#   <Group> = Chordata  Plantae  Arthropoda  Fungi  Bacteria  Archaea
#   (the NCBI taxon "Viridiplantae" is stored in the folder "Plantae")
#
# Usage (from anywhere):
#   bash Download_Data/scripts/download/download_cellular.sh                 # all groups
#   bash Download_Data/scripts/download/download_cellular.sh Archaea Fungi   # only these
#
# Optional environment variables:
#   CSPR_FASTA_GFF  output root   (default: <repo>/Fasta_GFF_inputs, same as config.py)
#   WORKERS         parallel downloads in the rehydrate step (default: 8)
#
# Safe to re-run: finished groups are skipped, interrupted groups are resumed
# (rehydrate only fetches the missing files).

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
OUT_ROOT="${CSPR_FASTA_GFF:-$REPO_ROOT/Fasta_GFF_inputs}"
WORKERS="${WORKERS:-8}"

DEFAULT_TAXA=(Chordata Viridiplantae Arthropoda Fungi Bacteria Archaea)
if [[ $# -gt 0 ]]; then TAXA=("$@"); else TAXA=("${DEFAULT_TAXA[@]}"); fi

for cmd in datasets unzip; do
  command -v "$cmd" >/dev/null 2>&1 || {
    echo "ERROR: '$cmd' not found. See Download_Data/README.md (section 'Install NCBI Datasets on Ubuntu')." >&2
    exit 1
  }
done

mkdir -p "$OUT_ROOT/_dehydrated_zips"
cd "$OUT_ROOT" || exit 1

LOG="$OUT_ROOT/download_log.tsv"
[[ -f "$LOG" ]] || printf "date\ttaxon\tfolder\tdatasets_version\tassemblies\n" > "$LOG"
DATASETS_VERSION="$(datasets --version 2>/dev/null | head -n1)"

download_one() {
  local taxon="$1" folder="$1"
  [[ "$taxon" == "Viridiplantae" ]] && folder="Plantae"      # folder name used by the analysis scripts
  local zip="_dehydrated_zips/${folder}.zip"

  if [[ -f "$folder/.rehydrated.ok" ]]; then
    echo "[$taxon] already complete ($folder), skipping"
    return 0
  fi

  # Step 1: dehydrated package (metadata + file list only, small)
  if [[ ! -s "$zip" ]]; then
    echo "[$taxon] Step 1/3: downloading the dehydrated package"
    datasets download genome taxon "$taxon" --reference --annotated --assembly-level complete,chromosome --include genome,gff3 --dehydrated --filename "$zip" || return 1
  fi

  # Step 2: unzip into Fasta_GFF_inputs/<folder>/
  if [[ ! -d "$folder/ncbi_dataset" ]]; then
    echo "[$taxon] Step 2/3: unzipping"
    unzip -q -o "$zip" -d "$folder" || return 1
  fi

  # Step 3: rehydrate (downloads the real FASTA/GFF files; resumable)
  echo "[$taxon] Step 3/3: rehydrating with $WORKERS workers"
  datasets rehydrate --directory "$folder" --max-workers "$WORKERS" || return 1

  touch "$folder/.rehydrated.ok"
  local n
  n=$(find "$folder/ncbi_dataset/data" -mindepth 1 -maxdepth 1 -type d -name 'GC*' 2>/dev/null | wc -l)
  printf "%s\t%s\t%s\t%s\t%s\n" "$(date +%F)" "$taxon" "$folder" "$DATASETS_VERSION" "$n" >> "$LOG"
  echo "[$taxon] done: $n assemblies in $OUT_ROOT/$folder"
}

FAILED=()
for taxon in "${TAXA[@]}"; do
  if ! download_one "$taxon"; then
    echo "[$taxon] FAILED, continuing with the next group" >&2
    FAILED+=("$taxon")
  fi
done

echo
if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo "Failed groups: ${FAILED[*]}"
  echo "Re-run the same command; finished groups are skipped and failed ones resume."
  exit 1
fi
echo "All requested groups finished. Log: $LOG"
