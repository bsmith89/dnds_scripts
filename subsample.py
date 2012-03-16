#!/usr/bin/env python
"""Take fasta-form nucleotide seqs to stdin, prints subsample to stdout.

"""

from Bio import SeqIO
import sys
import random

seq_records = list(SeqIO.parse(sys.stdin, "fasta"))
n = int(sys.argv[1])
sub_sample = random.sample(seq_records, n)
SeqIO.write(sub_sample, sys.stdout, "fasta")
