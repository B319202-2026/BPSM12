#!/usr/bin/env python3
"""
BPSM12 Lecture 12 - Split genomic DNA into coding/non-coding
Author: s2943447
"""

import os
import shutil
from Bio import Entrez, SeqIO

# === Settings ===
Entrez.email = "s2943447@ed.ac.uk"
local_data = "/localdisk/data/BPSM/Lecture12/plain_genomic_seq.txt"
work_dir = os.path.expanduser("~/Exercises/Lecture12")
os.chdir(work_dir)

# === Helper: write FASTA ===
def write_fasta(header, seq, filename):
    seq = seq.upper()
    with open(filename, "w") as f:
        f.write(">" + header + " length=" + str(len(seq)) + "bp\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")
    print("  Wrote " + filename + " (" + str(len(seq)) + " bp)")

# ============================================================
# Local sequence: plain_genomic_seq.txt
# Exon 1 = bases 1-63, Exon 2 = bases 91-end, Intron = bases 64-90
# ============================================================
print("=" * 50)
print("Local sequence")
print("=" * 50)

# Copy the file
shutil.copy(local_data, "local_genomic.txt")
print("Copied " + local_data + " -> local_genomic.txt")

# Read sequence
with open("local_genomic.txt") as f:
    genomic = f.read().strip().upper()

print("Full length: " + str(len(genomic)) + " bp")

# Split (Python is 0-based)
exon1 = genomic[0:63]    # bases 1-63
intron = genomic[63:90]   # bases 64-90
exon2 = genomic[90:]      # bases 91 to end

coding = exon1 + exon2
noncoding = intron

print("  Exon 1: " + str(len(exon1)) + " bp")
print("  Intron: " + str(len(intron)) + " bp")
print("  Exon 2: " + str(len(exon2)) + " bp")
print("  Total coding: " + str(len(coding)) + " bp")
print("  Total non-coding: " + str(len(noncoding)) + " bp")

write_fasta("local_genomic_coding", coding, "local_genomic_coding.fasta")
write_fasta("local_genomic_noncoding", noncoding, "local_genomic_noncoding.fasta")

# ============================================================
# Remote sequence: AJ223353
# ============================================================
print()
print("=" * 50)
print("Remote sequence AJ223353")
print("=" * 50)

# Fetch GenBank record from NCBI (includes annotation)
handle = Entrez.efetch(db="nucleotide", id="AJ223353",
                       rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

full_seq = str(record.seq).upper()
print("Full length: " + str(len(full_seq)) + " bp")

# Find CDS features
cds_ranges = []
for feature in record.features:
    if feature.type == "CDS":
        start = int(feature.location.start) + 1  # convert to 1-based
        end = int(feature.location.end)
        cds_ranges.append((start, end))
        print("  CDS: bases " + str(start) + "-" + str(end))

# Mark which positions are coding
coding_mask = [False] * len(full_seq)
for cds_start, cds_end in cds_ranges:
    for i in range(cds_start - 1, cds_end):
        coding_mask[i] = True

# Build coding and non-coding strings
coding_parts = []
noncoding_parts = []
for i, base in enumerate(full_seq):
    if coding_mask[i]:
        coding_parts.append(base)
    else:
        noncoding_parts.append(base)

coding_remote = "".join(coding_parts)
noncoding_remote = "".join(noncoding_parts)

print("  Coding: " + str(len(coding_remote)) + " bp")
print("  Non-coding: " + str(len(noncoding_remote)) + " bp")

write_fasta("AJ223353_coding", coding_remote, "AJ223353_coding.fasta")
write_fasta("AJ223353_noncoding", noncoding_remote, "AJ223353_noncoding.fasta")

# ============================================================
# Done
# ============================================================
print()
print("=" * 50)
print("Done! Generated files:")
print("=" * 50)
for f in sorted(os.listdir(".")):
    if f.endswith(".fasta"):
        print("  " + f)
