#!/usr/bin/env python
"""Removes all leaves but those given in names-file.


"""
import optparse
import sys
import Bio
from Bio import Phylo


def prune_all_but(tree, names):
    all_leaves = tree.get_terminals()
    prune_list = []
    for leaf in all_leaves:
        prune_list += [leaf.name]
    for name in names:
        prune_list.remove(name)
    for name in prune_list:
        tree.prune(name)


def main():
    usage = "usage: %prog [options] [tree] [names]"
    parser = optparse.OptionParser(usage = usage)
    (opts, args) = parser.parse_args()
    tree_file = None
    names_file = None
    try:
        names_path = sys.argv[2]
    except IndexError:
        tree_file = sys.stdin
        names_path = sys.argv[1]
    else:
        tree_file = open(sys.argv[1])
    finally:
        names_file = open(names_path)
    names = []
    for line in names_file:
        names += [line.strip()]
    tree = Bio.Phylo.read(tree_file, 'newick')
    prune_all_but(tree, names)
    Bio.Phylo.write(tree, sys.stdout, 'newick')

if __name__ == "__main__":
    main()
