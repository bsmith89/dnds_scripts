#!/usr/bin/env python
"""Unroot a newick formatted tree.

Take any newick formatted tree to stdin and write an unrooted tree
to stdout.

"""
from Bio import Phylo
import sys

def main():
    tree = Phylo.read(sys.stdin, 'newick')
    unroot(tree)
    Phylo.write(tree, sys.stdout, 'newick')

def unroot(tree):
    assert len(tree.clade.clades) == 2
    clade0 = tree.clade.clades[0]
    clade1 = tree.clade.clades[1]
    subclade0 = clade1.clades[0]
    subclade1 = clade1.clades[1]
    clade0.branch_length = float(clade0.branch_length) +\
                           float(clade1.branch_length)
    tree.clade.clades = [clade0, subclade0, subclade1]

if __name__ == '__main__':
    main()
