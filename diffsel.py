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
    analysis = DSAnalysis(fn_path, design_path, **opts)
    analysis.run()
    
    
class DSAnalysis():
    """
    
    """
    rep_dir_template = "rep%02d"
    hyp_dir_template = "h%d"
    
    def __init__(self, fn_path, design_path, **kwargs):
        # parse everything in kwargs
        self.num_reps = kwargs['num_reps']
        self.num_seqs = kwargs['num_seqs']
        if self.num_reps is None:
            if self.num_seqs is None:
                self.num_reps = 1
        else:
            # if num_reps is not None, make sure that a number
            # of sequences per sample has been specified.
            assert self.num_seqs is not None
        self.name = kwargs['name']
        if self.name is None:
            self.name = str(datetime.now())
        self.analysis_dir = os.path.join(kwargs['top_dir'], self.name)
        self._init_file_tree()
    
    def _init_file_tree(self):
        for rep in range(self.num_reps):
            for hyp in [0,1]:
                directory = os.path.join(self.analysis_dir,
                                         self.rep_dir_template % rep,
                                         self.hyp_dir_template % hyp)
                os.makedirs(directory)
    
    def run(self):
        self._translate()
        self._align_codons()
        self._generate_tree()
        self._subsample()
        self._label_trees()
        self._convert_afns()
        self._run_paml
    
    def __str__(self):
        pass
    
    def __repr__(self):
        pass
    
    
        
if __name__ == "__main__":
    main()























