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
        sys.stderr.write("Found leaf %s\n" % clade.name)
        try:
            if design[clade.name] is target:
                label_leaf(clade, labels[target])
                return True
            else:
                return False
        # if the node is not in the design mapping, then it's
        # either been excluded or it's already been labeled.
        # Either way it isn't the target grouping.
        except KeyError:
            return False
    else: # if it's a node
        sys.stderr.write("Entered node %s\n" % repr(clade))
        all_in_group = True
        for subclade in clade.clades:
            if not recurse_label(subclade, target, design, labels):
                all_in_group = False
        if all_in_group is False:
            return False
        else:
            label_node(clade, labels[target])
            return True

def label_leaf(leaf, label):
    leaf.name += label

def label_node(node, label):
    node.confidence = label

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
    sys.stderr.write("%s\n" % str(design))
    if labels_file is not None:
        for mapping in labels_file:
            grouping, label = mapping.split()
            labels[grouping] = label
    sys.stderr.write("%s\n" % str(labels))
    # Now actually do the labeling...
    for grouping in labels:
        sys.stderr.write("Applying labels to group %s\n" % grouping)
        recurse_label(tree.clade, grouping, design, labels)
    sys.stderr.write("%s\n" % str(tree))
    Bio.Phylo.write(tree, sys.stdout, 'newick', format_confidence = '%s')

if __name__ == '__main__':
    main()
