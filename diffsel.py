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
import random
import threading
import Bio.Phylo
import Bio.SeqIO
import Bio.AlignIO
import translate
import backalign
import treeassign
import paretree
from Bio.Align.Applications import MafftCommandline
from Bio.Application import AbstractCommandline, _Switch, _Argument
from Bio.Phylo.PAML import codeml


MAFFT_EXE = "/home/bjsmith/bin/mafft"
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
                      default = 1,
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
#    parser.add_option("-D", "--dry-run", "--do-not-run-paml", action = "store_false",
#                      dest = "run_paml", default = True,
#                      help = "dry run.  Don't run paml, only setup all \
#                      of the necessary files.")
    (opts, args) = parser.parse_args()
    fn_path = args[0]
    design_path = args[1]
    afa_path = opts.afa_path
    tree_path = opts.tree_path
    out_path = opts.out_path
    out_file = None
    if out_path is None:
        out_file = sys.stdout
    else:
        out_file = open(out_file, 'w')
    top_dir = os.path.abspath(opts.top_dir)
    num_reps = opts.num_reps
    num_seqs = opts.num_seqs
    if num_seqs is None and num_reps != 1:
        raise ValueError("Num_seqs is None, therefore num_reps must be 1.  \
Instead, num_reps = %d" % num_reps)
    num_trees = 1
    if tree_path is not None:
        num_trees = len(list(Bio.Phylo.parse(tree_path)))
    name = opts.name
    if name is None:
        d = datetime.datetime.now()
        name = "%d-%d-%d-%d-%d-%d" % (d.year, d.month, d.day, d.hour, d.minute, d.second)
    analysis_dir = os.path.join(top_dir, name)
    init_filesys(analysis_dir, num_reps, num_trees)
    if afa_path is None:
        fa_path = os.path.join(analysis_dir, "master.fa")
        _translate(fn_path, fa_path)
        afa_path = os.path.join(analysis_dir, "master.afa")
        mafft_align(fa_path, afa_path)
    # make a FastTree tree or split the tree(s) found in the tree file
    if tree_path is None:
        tree_path = os.path.join(analysis_dir, "master.tre01.tre")
        fasttree(afa_path, tree_path)
    else:
        split_trees(tree_path, analysis_dir)
    afn_path = os.path.join(analysis_dir, "master.afn")
    _backalign(fn_path, afa_path, afn_path)
    new_design_path = os.path.join(analysis_dir, "master.design")
    shutil.copy2(design_path, new_design_path)
    design_path = new_design_path
    _subsample(design_path, afn_path, num_reps, num_trees, num_seqs, analysis_dir)
    make_labels(analysis_dir, design_path)
    label_trees(analysis_dir, design_path, num_reps, num_trees)
    fa2phy_all(analysis_dir, num_reps)
    thread_list = paml_run_all(analysis_dir, num_reps, num_trees)
    
def init_filesys(analysis_dir, num_reps, num_trees):
    """Create directory tree required for analysis.
    
    './[analysis_dir]/rep##/tre##/h#/' for each rep, tree, and h0/h1. 
    
    """
    try:
        os.makedirs(analysis_dir)
    except OSError:
        raise ValueError("The directory '%s' appears to have already been \
taken, ostensibly by a previous analysis.  Either change the name or delete \
the previous analysis directory; CAUTION: make sure the previous analysis \
has finished running or files could be corrupted." % (analysis_dir))
    for rep in range(num_reps):
        for tree in range(num_trees):
            for hypo in [0,1]:
                directory = os.path.join(analysis_dir, "rep%02d" % (rep + 1), 
                                         "tree%02d" % (tree + 1), 
                                         "h%d" % hypo)
                os.makedirs(directory)        

def split_trees(tree_path, analysis_dir):
    """Separate each of the trees in tree_path into their own file.
    
    Outputs each tree to '[analysis_dir]/master.tre##.tre'.
    """
    i = 0
    for tree in Bio.Phylo.parse(tree_path, 'newick'):
        tree_i_path = os.path.join(analysis_dir, "master.tre%02d.tre" % (i + 1))
        Bio.Phylo.write(tree, tree_i_path, 'newick')
        i += 1

def _translate(fn_path, fa_path):
    """Translate the fn_sequences.
    
    Takes unaligned nucleotide FASTA from fn_path and writes amino-acid
    FASTA to out_path.
    
    """
    translate.files(fn_path, fa_path)

def mafft_align(fa_path, afa_path):
    """Align amino acid FASTA file.
    
    Takes amino-acid seqs from fa_path and writes aligned amino-acids
    to afa_path. 
    
    """
    mafft_call = MafftCommandline(input = fa_path)
    mafft_call.maxiterate = 1000
    mafft_call.retree = 2
    stdout, stderr = mafft_call()
    open(afa_path, "w").write(stdout)
    open("%s.err" % afa_path, 'w').write(stderr)

def _backalign(fn_path, afa_path, out_path):
    """Backalign the fn sequences to match the afa alignment.
    
    Takes nucleotide FASTA file from fn_path, aligned amino-acid FASTA
    at afa_path, and prints aligned nucleotide FASTA to out_path.
    
    """
    backalign.files(fn_path, afa_path, out_path)

def fasttree(afa_path, tree_path):
    """Construct a phylogenetic tree using FastTree.
    
    Takes an amino-acid alignment from afa_path and outputs a
    newick-formatted, support-less tree using FastTree to tree_path.
    
    """
    class FastTreeCommandline(AbstractCommandline):
        def __init__(self, cmd="FastTreeMP", **kwargs):
            self.parameters = []
            self.parameters += \
                [_Switch(['-nosupport', 'nosupport'],
                         "don't include support values in the output tree")]
            self.parameters += \
                [_Argument(['', 'input'],
                           'input file')]
            # set a list of parameters which are objects derived from the base class
            # _AbstractParameter
            AbstractCommandline.__init__(self, cmd, **kwargs)
    fasttree_call = FastTreeCommandline(input = afa_path)
    fasttree_call.nosupport = True
    stdout, stderr = fasttree_call()
    open(tree_path, 'w').write(stdout)
    open("%s.err" % tree_path, 'w').write(stderr)

        
def _subsample(design_path, afn_path, num_reps, num_trees, num_seqs, analysis_dir):
    """Subsample num_seqs seqs from each treatment in design_path.
    
    (1) Subsamples num_seqs, num_reps times from the names presented in
    design_path and moves these aligned nucleotide sequences to 
    '[analysis_dir]/rep##/rep##.afn'.
    (2) Pares down the tree found in each master.tre##.tre based
    on the sequences in the subsample, and outputs the result to
    '[analysis_dir]/rep##/rep##.afn'.
    
    """
    for i in range(num_reps):
        afn_out_path = os.path.join(analysis_dir,
                                    "rep%02d" % (i + 1),
                                    "rep%02d.afn" % (i + 1))
        single_subsample(design_path, afn_path, num_seqs, afn_out_path)
        for j in range(num_trees):
            master_tree_path = os.path.join(analysis_dir,
                                            "master.tre%02d.tre" % 
                                            (j + 1))
            out_tree_path = os.path.join(analysis_dir, "rep%02d" % (i + 1),
                                         "tree%02d" % (j + 1),
                                         "rep%02d.tre%02d.tre" % (i + 1, j + 1))
            names_list = names_from_fa(afn_out_path)
            _paretree(master_tree_path, names_list, out_tree_path)

def names_from_fa(fa_path):
    names = []
    seq_records = Bio.SeqIO.parse(fa_path, "fasta")
    for record in seq_records:
        names += [record.name]
    return names
        
def _paretree(master_tree_path, names_list, out_tree_path):
    paretree.files(master_tree_path, names_list, out_tree_path)

def single_subsample(design_path, afn_path, num_seqs, afn_out_path):
    """Make a single subsample of seqs in afn_path.
    
    Subsamples from afn_path and outputs the resulting sequences to
    afn_out_path.
    
    """
    if num_seqs is None:
        shutil.copy2(afn_path, afn_out_path)
        return
    
    full_design = {}
    with open(design_path) as design_file:
        for line in design_file:
            name, treatment = line.strip().split()
            full_design[name] = treatment
    treatments = {}
    for treatment in set(full_design.values()):
        treatments[treatment] = []
    seq_records = list(Bio.SeqIO.parse(afn_path, 'fasta'))
    for record in seq_records:
        treatments[full_design[record.name]] += [record.name]
    sampled_seq_names = []
    for treatment in treatments:
        sampled_seq_names += random.sample(treatments[treatment], num_seqs)
    out_seqs = []
    for record in seq_records:
        if record.name in sampled_seq_names:
            out_seqs += [record]
    Bio.SeqIO.write(out_seqs, afn_out_path, 'fasta')

def make_labels(analysis_dir, design_path):
    """Take a design file and make appropriate h0 and h1 labels.
    
    The labels are output to '[analysis_dir]/master.h#.labels'.
    """
    h0_labels_path = os.path.join(analysis_dir, "master.h0.labels")
    h1_labels_path = os.path.join(analysis_dir, "master.h1.labels")
    h0_labels_file = open(h0_labels_path, 'w')
    h1_labels_file = open(h1_labels_path, 'w')
    design_file = open(design_path)
    treatments = []
    for line in design_file:
        treatment = line.strip().split()[1]
        if treatment not in treatments:
            treatments += [treatment]
    design_file.close()
    i = 0
    for treatment in treatments:
        h0_labels_file.write("%s\t#1\n" % (treatment))
        h1_labels_file.write("%s\t#%d\n" % (treatment, i + 1))
        i += 1
    h0_labels_file.close()
    h1_labels_file.close()

def label_trees(analysis_dir, design_path, num_reps, num_trees):
    """Make PAML labeled trees in all of the appropriate directories.
    
    Takes master trees from '[analysis_dir]/master.tre##.tre' and
    outputs '[analysis_dir]/rep##/tre##/h#/rep##.tre##.h#.nwk'
    based on the labels found in '[analysis_dir]/master.h#.labels'.
    
    """
    for rep in range(num_reps):
        for tree in range(num_trees):
            for hypo in [0,1]:
                out_path = os.path.join(analysis_dir, "rep%02d" % (rep + 1), 
                                        "tree%02d" % (tree + 1), 
                                        "h%d" % hypo,
                                        "rep%02d.tre%02d.h%d.nwk" % \
                                        (rep + 1, tree + 1, hypo))
                tree_path = os.path.join(analysis_dir, "rep%02d" % (rep + 1),
                                         "tree%02d" % (tree + 1),
                                         "rep%02d.tre%02d.tre" % \
                                         (rep + 1, tree + 1))
                label_path = os.path.join(analysis_dir, 
                                          "master.h%d.labels" % \
                                          (hypo))
                treeassign.files(tree_path, design_path, label_path, out_path)

def fa2phy_all(analysis_dir, num_reps):
    """Convert the aligned FASTA files to Phylip format.
    
    Makes them in the appropriate PAML format (interleaved,
    with the 'I' option) and outputs to
    '[analysis_dir]/rep##/rep##.phy'
    """
    for rep in range(num_reps):
        fa_path = os.path.join(analysis_dir, "rep%02d" % (rep + 1),
                               "rep%02d.afn" % (rep + 1))
        fa_file = open(fa_path)
        phy_path = os.path.join(analysis_dir, "rep%02d" % (rep + 1),
                               "rep%02d.phy" % (rep + 1))
        phy_file = open(phy_path, 'w')
        alignments = Bio.AlignIO.parse(fa_file, 'fasta')
        Bio.AlignIO.write(alignments, phy_file, 'phylip')
        phy_file.close()
        phy_file = open(phy_path)
        phy_lines = phy_file.readlines()
        phy_file.close()
        phy_file = open(phy_path, 'w')
        first = True
        for line in phy_lines:
            if first:
                phy_file.write("%s I\n" % line.strip())
                first = False
            else:
                phy_file.write(line)
        phy_file.close()
        
def paml_run_all(analysis_dir, num_reps, num_trees):
    thread_list = []
    for i in range(num_reps):
        for j in range(num_trees):
            for k in [0,1]:
                cml = codeml.Codeml()
                cml.working_dir = os.path.join(analysis_dir,
                                               "rep%02d" % (i + 1),
                                               "tree%02d" % (j + 1),
                                               "h%d" % k)
                cml.alignment = os.path.join(analysis_dir,
                                             "rep%02d" % (i + 1),
                                             "rep%02d.phy" % (j + 1))
                cml.tree = os.path.join(cml.working_dir,
                                        "rep%02d.tre%02d.h%d.nwk" % (i + 1, j + 1, k))
                cml.ctl_file = os.path.join("rep%02d.tre%02d.h%d.ctl" % (i + 1, j + 1, k))
                cml.out_file = os.path.join(analysis_dir,
                                        "rep%02d" % (i + 1),
                                        "tree%02d" % (j + 1),
                                        "h%d" % k,
                                        "rep%02d.tre%02d.h%d.mlc" % (i + 1, j + 1, k))
                cml.set_options(runmode = 0,
                                seqtype = 1,
                                CodonFreq = 2,
                                ndata = 1,
                                clock = 0,
                                aaDist = 0,
                                model = 2,
                                NSsites = [0],
                                icode = 0,
                                Mgene = 0,
                                fix_kappa = 0,
                                kappa = 3.53330,
                                fix_omega = 0,
                                omega = 1,
                                getSE = 0,
                                RateAncestor = 0,
                                Small_Diff = 0.0000005,
                                cleandata = 0,
                                fix_blength = 1,
                                method = 0)
                cml.write_ctl_file()
                thread = CodemlThread(cml)
                print("Adding a thread")
                thread.start()
                thread_list.append(thread)
    return thread_list
                
class CodemlThread(threading.Thread):
    def __init__(self, cml_object):
        self.cml = cml_object
        threading.Thread.__init__(self)
        
    def run(self):
        print(os.path.abspath(self.cml.tree))
        self.cml.run()
    
                        
if __name__ == "__main__":
    main()
