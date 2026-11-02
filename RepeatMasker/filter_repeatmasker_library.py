"""
filter_repeatmasker_library.py -- list the distinct taxon names that occur in a Dfam FASTA
library (second word of each header), to help choose which lineage-specific entries to keep
when building a custom RepeatMasker library.  Prints only; changes nothing.

Usage:
    python3 filter_repeatmasker_library.py /path/to/dfam-fasta-download.fasta
"""
import argparse
from Bio import SeqIO

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("library", help="Dfam FASTA library file")
a = ap.parse_args()

seen = set()
for te in SeqIO.parse(a.library, "fasta"):
    parts = te.description.split()
    if len(parts) < 2 or parts[1] in seen:
        continue
    seen.add(parts[1])
    print(len(seen), parts[1])
