#!/usr/bin/env python
"""Add identifiers to single-treatment clades.

Takes a newick formatted tree to stdin, a design file argument
and an optional labels file argument, and prints a labeled tree.

Every leaf which belongs to a particular treatment gets a tag,
as do all higher level nodes which have no sub-clades containing
leaves which belong to another treatment group.  Leaves which are
not assigned to a treatment are ignored, and do not count as either
in or out-of-treatment clades.

"""

import Bio
from Bio import Phylo
import sys

def recurse_label_nodes(clade, target, design, labels):
    """Label nodes and return True if all children are of Group.

    """
    if clade.is_terminal(): # if it's a leaf
        try:
            if design[clade.name] == target:
                return True
            else:
                return False
        # If the node is not in the design mapping ignore it
        # and don't consider it in marking higher level nodes.
        except KeyError:
            return None
    else: # if it's a node
        all_in_group = None
        for subclade in clade.clades:
            subclade_in_group = recurse_label_nodes(subclade, target, design, labels)
            if subclade_in_group is False:
                all_in_group = False
            elif subclade_in_group is True and all_in_group is not False:
                all_in_group = True
        if all_in_group is False:
            return False
        elif all_in_group is True:
            label_node(clade, labels[target])
            return True
        else:
            return None

def label_all(tree, design, labels):
    for grouping in labels:
        recurse_label_nodes(tree.clade, grouping, design, labels)
        # Because labeling directly changes the name, I can't label
        # leaves before I'm done with all of the nodes.
    for leaf in tree.get_terminals():
        try:
            label_leaf(leaf, labels[design[leaf.name]])
        except KeyError:
            continue

def label_leaf(leaf, label):
    """Label the given leaf with the provided label string.
    
    Adds the label, after a space, to the leaf's name.
    """
    leaf.name += " %s" % label

def label_node(node, label):
    """Label the given node with the provided label string.

    Adds the label, with a space, as the confidence value on the tree.
    This is a hack, since Bio.Phylo doesn't seem to have another way to
    label a node or to print this label the way I want it to a newick
    tree.
    
    To correctly print a tree with these node labels you'll need to use
    Bio.Phylo.write(tree, file, format, *format_confidence = '%s'*) since
    otherwise it outputs confidence as decimal padded float.

    """
    node.confidence = " %s" % label

def parse_design_and_labels(design_file, labels_file = None):
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
    return design, labels


def write_tree(tree, out_file):
    Bio.Phylo.write(tree, out_file, 'newick', format_confidence = '%s')

def main():
    design_path = sys.argv[1]
    design_file = open(design_path)
    try:
        labels_path = sys.argv[2]
        labels_file = open(labels_path)
    except IndexError:
        labels_file = None
    tree = Bio.Phylo.read(sys.stdin, 'newick')
    design, labels = parse_design_and_labels(design_file, labels_file)
    label_all(tree, design, labels)
    write_tree(tree, sys.stdout)

if __name__ == '__main__':
    main()
