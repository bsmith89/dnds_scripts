#!/usr/bin/env python

import Bio
from Bio import SeqIO
import sys

seq_records = list(Bio.SeqIO.parse(sys.stdin, 'fasta'))
Bio.SeqIO.write(seq_records, sys.stdout, 'phylip')
