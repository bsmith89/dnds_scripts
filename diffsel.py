#!/usr/bin/env python
"""Test differential dN/dS for lineages unique to treatments.

A pipelining of my dN/dS analysis.  Takes a design file, nucleotide
sequence, and optional amino-acid alignment in FASTA format and tests
for differential dN/dS between lineages present in one treatment vs.
the other.

TODO: Consider re-writing this as a analysis object which stores all
of the important variables.

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
    fa2phy_all(analysis_dir, num_reps)
    paml_run_all(analysis_dir)
    
    
def paml_run_all(analysis_dir, out_file):
    """Carry out the full PAML analysis and print result to out_file.
    
    """
    pass

def fa2phy_all(analysis_dir, num_reps):
    """Convert the aligned FASTA files to Phylip format.
    
    Makes them in the appropriate PAML format (interleaved,
    with the 'I' option) and outputs to
    '[analysis_dir]/rep##/name.rep##.phy'
    """
    pass
    
def make_labels(analysis_dir, name, design_path):
    """Take a design file and make appropriate h0 and h1 labels.
    
    The labels are output to '[analysis_dir]/[name].master.h#.labels'.
    """
    pass
                
                
def label_trees(analysis_dir, design_path, num_reps, num_trees):
    """Make PAML labeled trees in all of the appropriate directories.
    
    Takes master trees from '[analysis_dir]/[name].master.tre##.tre' and
    outputs '[analysis_dir]/rep##/tre##/h#/[name].rep##.tre##.h#.nwk'
    based on the labels found in '[analysis_dir]/[name].master.h#.labels'.
    
    """
    for rep in range(num_reps):
        for tree in range(num_trees):
            for hypothesis in [0,1]:
                directory = os.path.join(analysis_dir, "rep%02d" % (rep + 1), 
                                         "tree%02d" % (tree + 1), 
                                         "h%d" % hypothesis)
                pass
        
def split_trees(tree_path, analysis_dir):
    """Separate each of the trees in tree_path into their own file.
    
    Outputs each tree to '[analysis_dir]/[name].master.tre##.tre'.
    """
    pass
        
def subsample(design_path, afn_path, num_reps, num_trees, num_seqs, analysis_dir):
    """Subsample num_seqs seqs from each treatment in design_path.
    
    (1) Subsamples num_seqs, num_reps times from the names presented in
    design_path and moves these aligned nucleotide sequences to 
    '[analysis_dir]/rep##/name.rep##.afn'.
    (2) Pares down the tree found in each [name].master.tre##.tre based
    on the sequences in the subsample, and outputs the result to
    '[analysis_dir]/rep##/name.rep##.afn'.
    
    """
    pass

def backalign(fn_path, afa_path, out_path):
    """Backalign the fn sequences to match the afa alignment.
    
    Takes nucleotide FASTA file from fn_path, aligned amino-acid FASTA
    at afa_path, and prints aligned nucleotide FASTA to out_path.
    
    """
    pass

def translate(fn_path, fa_path):
    """Translate the fn_sequences.
    
    Takes unaligned nucleotide FASTA from fn_path and writes amino-acid
    FASTA to out_path.
    
    """
    pass

def mafft_align(fa_path, afa_path):
    """Align amino acid FASTA file.
    
    Takes amino-acid seqs from fa_path and writes aligned amino-acids
    to afa_path. 
    
    """
    pass

def fasttree(afa_path, tree_path):
    """Construct a phylogenetic tree using FastTree.
    
    Takes an amino-acid alignment from afa_path and outputs a
    newick-formatted, support-less tree using FastTree to tree_path.
    
    """
    pass

def init_filesys(analysis_dir, num_reps, num_trees):
    """Create directory tree required for analysis.
    
    './[analysis_dir]/rep##/tre##/h#/' for each rep, tree, and h0/h1. 
    
    """
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
