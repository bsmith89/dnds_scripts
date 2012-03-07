#!/usr/bin/env python
"""

"""

import Bio
from Bio import Phylo
import sys

def recurse_label(clade, target, design, labels):
    """Label nodes and return True if all children are of Group.

    Just pseudo-code right now
    """
    if clade.is_terminal(): # if it's a leaf
        try:
            if design[clade.name] is target:
                return True
            else:
                return False
        # If the node is not in the design mapping ignore it
        # and don't consider it in marking higher level nodes.
        except KeyError:
            return None # or should it be True? nothing at all?
    else: # if it's a node
        all_in_group = True
        for subclade in clade.clades:
            if recurse_label(subclade, target, design, labels) is False:
                all_in_group = False
        if all_in_group is False:
            return False
        else:
            label_node(clade, labels[target])
            return True

def label_leaf(leaf, label):
    leaf.name += " %s" % label

def label_node(node, label):
    node.confidence = " %s" % label

def main():
    design_path = sys.argv[1]
    design_file = open(design_path)
    try:
        labels_path = sys.argv[2]
        labels_file = open(labels_path)
    except IndexError:
        labels_file = None
    tree = Bio.Phylo.read(sys.stdin, 'newick')
    design = {}
    labels = {}
    for mapping in design_file:
        item, grouping = mapping.split()
        design[item] = grouping
        if labels_file is None and grouping not in labels:
            labels[grouping] = "#%d" % (len(labels) + 1)
    if labels_file is not None:
        for mapping in labels_file:
            grouping, label = mapping.split()
            labels[grouping] = label
    # Now actually do the labeling...
    for grouping in labels:
        recurse_label(tree.clade, grouping, design, labels)
    # Because labeling directly changes the name, I can't label
    # leaves before I'm done with all of the nodes.
    for leaf in tree.get_terminals():
        try:
            label_leaf(leaf, labels[design[leaf.name]])
        except KeyError:
            continue
    Bio.Phylo.write(tree, sys.stdout, 'newick', format_confidence = '%s')

if __name__ == '__main__':
    main()
