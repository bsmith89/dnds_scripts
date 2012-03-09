#!/usr/bin/env python
"""Take fasta-form nucleotide seqs to stdin, print AA-seqs to stdout.

"""

from Bio import SeqIO
import sys

seq_records = list(SeqIO.parse(sys.stdin, "fasta"))
for record in seq_records:
    seq = record.seq
    aa = seq.translate()
    record.seq = aa
SeqIO.write(seq_records, sys.stdout, "fasta")
