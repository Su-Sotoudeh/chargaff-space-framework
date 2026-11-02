#!/usr/bin/env bash
# Report what was downloaded and flag problems.
#   bash Download_Data/scripts/download/check_downloads.sh

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
IN_ROOT="${CSPR_FASTA_GFF:-$REPO_ROOT/Fasta_GFF_inputs}"

echo "=== NCBI Datasets packages (assemblies, with FASTA, with GFF) ==="
shopt -s nullglob
for d in "$IN_ROOT"/*/ncbi_dataset/data; do
  group="$(basename "$(dirname "$(dirname "$d")")")"
  total=0; with_fna=0; with_gff=0; missing=()
  for a in "$d"/GC*/; do
    total=$((total+1))
    fna=("$a"*.fna); gff=("$a"genomic.gff)
    [[ ${#fna[@]} -gt 0 ]] && with_fna=$((with_fna+1))
    if [[ -s "${gff[0]:-}" ]]; then with_gff=$((with_gff+1)); else missing+=("$(basename "$a")"); fi
  done
  printf "%-28s total=%-6s fasta=%-6s gff=%-6s\n" "$group" "$total" "$with_fna" "$with_gff"
  if [[ ${#missing[@]} -gt 0 ]]; then
    printf '%s\n' "${missing[@]}" > "$IN_ROOT/${group}_missing_gff.txt"
    echo "   -> ${#missing[@]} assemblies without GFF listed in ${group}_missing_gff.txt"
  fi
done

echo
echo "=== RefSeq FTP files, optional (gzip integrity) ==="
for set in mitochondrion plastid plasmid viral; do
  dir="$IN_ROOT/RefSeq_FTP/$set"
  [[ -d "$dir" ]] || { echo "$set: not downloaded"; continue; }
  n=0; bad=0
  for f in "$dir"/*.gz; do
    n=$((n+1))
    gzip -t "$f" 2>/dev/null || { bad=$((bad+1)); echo "   corrupt: $f"; }
  done
  printf "%-14s files=%-5s corrupt=%s\n" "$set" "$n" "$bad"
done
