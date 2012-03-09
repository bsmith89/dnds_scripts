#!/usr/bin/env python
"""Test differential dN/dS for lineages unique to treatments.

A pipelining of my dN/dS analysis.  Takes a design file, newick
formatted tree, and nucleotide alignment and tests for any difference
between treatments as specified in the design file.

"""
import optparse
import sys
import os

def main():
    usage = "usage: %prog [options] tree aligned-nuc design"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option("-o", "--out", dest = "out_file",
                      help = "if given will print output here instead of stdout")
    parser.add_option("-d", "--dir", dest = "analysis_dir",
                      help = "if given will save intermediate files here instead of ./diffsel")
    (opts, args) = parser.parse_args()
    tree_file = open(args[0])
    afn_file = open(args[1])
    design_file = open(args[2])
    out_file = None
    if opts.out_file is not None:
        out_file = open(opts.out_file)
    else:
        out_file = sys.stdout
    adir = None
    if opts.analysis_dir is not None:
        adir = opts.analysis_dir
    else:
        adir = './diffsel/'
    # Now some pseudocode
#    create_dir(adir)
#    change_dir(adir)
#    write_ctl_file(**any_kwargs_from_options)
#    create_dir(h0_dir)
#    create_dir(h1_dir)
#    write_h0_labels()
    
    

if __name__ == "__main__":
    main()