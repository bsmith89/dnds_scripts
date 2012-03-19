#!/usr/bin/env python
"""Test differential dN/dS for lineages unique to treatments.

A pipelining of my dN/dS analysis.  Takes a design file, nucleotide
sequence, and optional amino-acid alignment in FASTA format and tests
for differential dN/dS between lineages present in one treatment vs.
the other.

"""
import optparse
import sys
import os
import shutil
import datetime
import Bio
from Bio import Phylo

DEFAULT_TOP_DIR = "./.tmpdiffsel"

def main():
    usage = "usage: %prog [options] fn design"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option("-a", "--afa", "--aligned-aa", dest = "afa_path",
                      default = None,
                      help = "amino-acid alignment in FASTA format.  \
                              DEFAULT: generates mafft alignment from \
                              translated nucleotides.")
    parser.add_option("-t", "--tree", dest = "tree_path", default = None,
                      help = "phylogenetic tree in newick format.  \
                              If multiple trees are found repeats all \
                              analysis with each tree; uses the same \
                              subsamples for each if reps > 1 \
                              DEFAULT: generates FastTree phylogeny from \
                              aligned amino-acids.")
    parser.add_option("-o", "--out", dest = "out_path", default = None,
                      help = "print output here instead of stdout")
    parser.add_option("-d", "--dir", dest = "top_dir",
                      default = DEFAULT_TOP_DIR,
                      help = "top level directory for intermediate analysis \
                              files.  DEFAULT: ./.tmpdiffsel")
    parser.add_option("-r", "--reps", dest = "num_reps", type = 'int',
                      DEFAULT = 1,
                      help = "number of replicate subsamplings.  \
                              DEFAULT: 1")
    parser.add_option("-n", "--subsample-size", dest = "num_seqs",
                      type = 'int', default = None,
                      help = "number of sequences in the subsample.  \
                              Increased sample size increases both \
                              power and analysis time.  DEFAULT: use \
                              all sequences.")
    parser.add_option("-N", "--name", "--analysis-name", dest = "name",
                      default = None,
                      help = "the name of the analysis directory under \
                              [top-dir]/.")
    (opts, args) = parser.parse_args()
    fn_path = args[0]
    design_path = args[1]
    afa_path = opts['afa_path']
    tree_path = opts['tree_path']
    out_path = opts['out_path']
    out_file = None
    if out_path is None:
        out_file = sys.stdout
    else:
        out_file = open(out_file, 'w')
    top_dir = opts['top_dir']
    num_reps = opts['num_reps']
    num_seqs = opts['num_seqs']
    if num_seqs is None and num_reps != 1:
        raise ValueError("Num_seqs is None, therefore num_reps must be 1.  \
                          Instead, num_reps = %d" % num_reps)
    num_trees = 1
    if tree_path is not None:
        num_trees = len(list(Bio.Phylo.parse(tree_path)))
    name = opts['name']
    if name is None:
        name = str(datetime.now())
    analysis_dir = os.path.join(top_dir, name)
    init_filesys(analysis_dir, num_reps, num_trees)
    if afa_path is None:
        fa_path = os.path.join(analysis_dir, "%s.master.fa" % name)
        translate(fn_path, fa_path)
        afa_path = os.path.join(analysis_dir, "%s.master.afa" % name)
        mafft_align(fa_path, afa_path)
    # make a FastTree tree or split the tree(s) found in the tree file
    if tree_path is None:
        tree_path = os.path.join(analysis_dir, "%s.master.tre01.tre" % name)
        fasttree(afa_path, tree_path)
    else:
        split_trees(tree_path, analysis_dir)
    afn_path = os.path.join(analysis_dir, "%s.master.afn" % name)
    backalign(fn_path, afa_path, afn_path)
    new_design_path = os.path.join(analysis_dir, "%s.master.design" % name)
    shutil.copy2(design_path, new_design_path)
    design_path = new_design_path
    subsample(design_path, afn_path, num_reps, num_trees, num_seqs, analysis_dir)
    make_labels(analysis_dir, design_path)
    label_trees(analysis_dir, design_path, num_reps, num_trees)
    
    
def make_labels(analysis_dir, design_path):
    pass
                
                
def label_trees(analysis_dir, design_path, num_reps, num_trees):
    for rep in range(num_reps):
        for tree in range(num_trees):
            for hypothesis in [0,1]:
                directory = os.path.join(analysis_dir, "rep%02d" % (rep + 1), 
                                         "tree%02d" % (tree + 1), 
                                         "h%d" % hypothesis)
                pass
        
def split_trees(tree_path, analysis_dir):
    pass
        
def subsample(design_path, afn_path, tree_path, num_reps, num_seqs, analysis_dir):
    pass

def backalign(fn_path, afa_path, out_path):
    pass

def translate(fn_path, fa_path):
    pass

def mafft_align(fa_path, afa_path):
    pass

def fasttree(afa_path, tree_path):
    pass

def init_filesys(analysis_dir, num_reps, num_trees):
    try:
        os.mkdirs(analysis_dir)
    except OSError:
        raise ValueError("The directory '%s' appears to have already been \
                          taken, ostensibly by a previous analysis.  Either \
                          change the name or delete the previous analysis \
                          directory; CAUTION: make sure the previous analysis \
                          has finished running or files could be corrupted." % 
                          (analysis_dir))
    for rep in range(num_reps):
        for tree in range(num_trees):
            for hypothesis in [0,1]:
                directory = os.path.join(analysis_dir, "rep%02d" % (num_reps + 1), 
                                         "tree%02d" % (num_trees + 1), 
                                         "h%d" % hypothesis)
                os.makedirs(directory)
        
if __name__ == "__main__":
    main()
