#!/usr/bin/env python
"""Converts fasta formatted stding to phylip-form stdout.

"""

import Bio
from Bio import SeqIO
import sys

seq_records = list(Bio.SeqIO.parse(sys.stdin, 'fasta'))
Bio.SeqIO.write(seq_records, sys.stdout, 'phylip')
