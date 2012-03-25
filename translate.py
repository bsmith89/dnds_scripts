#!/usr/bin/env python
"""Take fasta-form nucleotide seqs to stdin, print AA-seqs to stdout.

"""

from Bio import SeqIO
import sys

def files(in_path, out_path):
    seq_records = list(SeqIO.parse(in_path, "fasta"))
    _translate(seq_records)
    SeqIO.write(seq_records, out_path, "fasta")
    
def _translate(seq_records):
    for record in seq_records:
        seq = record.seq
        aa = seq.translate()
        record.seq = aa

def main():
    seq_records = list(SeqIO.parse(sys.stdin, "fasta"))
    for record in seq_records:
        seq = record.seq
        aa = seq.translate()
        record.seq = aa
    SeqIO.write(seq_records, sys.stdout, "fasta")

if __name__ == "__main__":
    main()
