#!/usr/bin/env bash
# Download RefSeq mitochondrion, plastid (chloroplast), plasmid and viral
# sequence files (FASTA + GenBank flat files) from the NCBI FTP server.
#
# The script reads the server's directory listing, so it does not depend on a
# hard-coded number of file parts.
#
# Usage:
#   bash Download_Data/scripts/download/download_organelles_plasmids_viruses.sh                  # all four sets
#   bash Download_Data/scripts/download/download_organelles_plasmids_viruses.sh mitochondrion    # one set
#   sets: mitochondrion  plastid  plasmid  viral
#
# Optional environment variable:
#   CSPR_FASTA_GFF  output root (default: <repo>/Fasta_GFF_inputs, same as config.py)
#
# Output: Fasta_GFF_inputs/RefSeq_FTP/<set>/*.genomic.{fna,gbff}.gz  (GenBank format).
# NOTE: GFF_Mapping/Mitochondrion and GFF_Mapping/Chloroplast read the NCBI Datasets package
# layout (Fasta_GFF_inputs/Mitochondrion|Chloroplast/ncbi_dataset/data/...), which is
# downloaded from the NCBI web page -- see Download_Data/README.md.  This script is an
# optional, reproducible alternative source (RefSeq release files).
#
# Safe to re-run: complete files are kept, partial files are resumed.

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
OUT_ROOT="${CSPR_FASTA_GFF:-$REPO_ROOT/Fasta_GFF_inputs}/RefSeq_FTP"
BASE="${NCBI_FTP_BASE:-https://ftp.ncbi.nlm.nih.gov}"   # override only for testing

declare -A URL=(
  [mitochondrion]="$BASE/refseq/release/mitochondrion/"
  [plastid]="$BASE/refseq/release/plastid/"
  [plasmid]="$BASE/genomes/refseq/plasmid/"
  [viral]="$BASE/refseq/release/viral/"
)
DEFAULT_SETS=(mitochondrion plastid plasmid viral)
if [[ $# -gt 0 ]]; then SETS=("$@"); else SETS=("${DEFAULT_SETS[@]}"); fi

for cmd in curl gzip; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: '$cmd' not found" >&2; exit 1; }
done

fetch_set() {
  local set="$1" url="${URL[$1]:-}"
  [[ -n "$url" ]] || { echo "Unknown set '$set' (use: ${DEFAULT_SETS[*]})" >&2; return 1; }

  local dir="$OUT_ROOT/$set"
  mkdir -p "$dir"

  echo "[$set] reading file list from $url"
  local listing files
  listing="$(curl -fsSL --retry 5 "$url")" || { echo "[$set] cannot read $url" >&2; return 1; }
  files="$(grep -oE 'href="[^"]+\.genomic\.(fna|gbff)\.gz"' <<<"$listing" | sed -E 's/^href="//; s/"$//' | sort -u)"
  [[ -n "$files" ]] || { echo "[$set] no *.genomic.fna.gz / *.genomic.gbff.gz files found; the server layout may have changed" >&2; return 1; }

  local total bad=0 f
  total=$(wc -l <<<"$files")
  echo "[$set] $total files to fetch into $dir"

  while IFS= read -r f; do
    if [[ -s "$dir/$f" ]] && gzip -t "$dir/$f" 2>/dev/null; then
      continue                                  # already complete and valid
    fi
    curl -fsSL --retry 5 --retry-delay 5 -C - -o "$dir/$f" "$url$f" || { echo "[$set] download failed: $f (re-run to resume)" >&2; bad=$((bad+1)); continue; }
    gzip -t "$dir/$f" 2>/dev/null || { echo "[$set] corrupt file removed: $f (re-run to fetch it again)" >&2; rm -f "$dir/$f"; bad=$((bad+1)); }
  done <<<"$files"

  # Provenance: what was downloaded, when, and checksums
  ( cd "$dir" && { echo "# $set downloaded $(date -u +%FT%TZ) from $url"; sha256sum ./*.gz; } > MANIFEST.txt )

  if [[ $bad -gt 0 ]]; then echo "[$set] finished with $bad problem file(s)" >&2; return 1; fi
  echo "[$set] done"
}

FAILED=()
for s in "${SETS[@]}"; do
  fetch_set "$s" || FAILED+=("$s")
done

echo
if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo "Sets with problems: ${FAILED[*]}  (re-run to resume)"
  exit 1
fi
echo "All requested sets finished. Output: $OUT_ROOT"
