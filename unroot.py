#!/usr/bin/env python
"""Unroot a newick formatted tree.

Take any newick formatted tree to stdin and write an unrooted tree
to stdout.

"""
from Bio import Phylo
import sys

tree = Phylo.read(sys.stdin, 'newick')
unrooted = tree.unroot()
Phylo.write(sys.stdout, 'newick')
