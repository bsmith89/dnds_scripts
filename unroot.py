#!/usr/bin/env python

from Bio import Phylo
import sys

tree = Phylo.read(sys.stdin, 'newick')
unrooted = tree.unroot()
Phylo.write(sys.stdout, 'newick')
