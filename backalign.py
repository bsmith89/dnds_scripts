#!/usr/bin/env python
"""Correctly aligns nucleotide sequences (codon-wise).

Takes AA-alignments with their associated nucleotide sequences and
outputs the aligned nucleotides.

TODO: Figure out how to deal with AA seqs that aren't the correct
length.
"""
import optparse
import sys
import Bio
from Bio import SeqIO, Alphabet, Seq

def get_seq_record(ident, seq_recs):
    """Returns the first seq record with the provided name."""
    for seq_rec in seq_recs:
        if seq_rec.id == ident:
            return seq_rec

def back_align_all(nucl_recs, aa_recs):
    out_nucl_recs = []
    for nucl_rec in nucl_recs:
        aa_rec = get_seq_record(ident = nucl_rec.id, seq_recs = aa_recs)
        try:
            nucl_rec.seq = back_align(nucl_rec, aa_rec)
        except ValueError as back_align_error:
            sys.stderr.write("WARNING: nucl_seq '%s' and aa_seq '%s' could not be back-aligned:\n%s\nSequences removed from analysis.\n" % (nucl_rec.id, aa_rec.id, str(back_align_error)))
            nucl_rec.seq = Bio.Seq.Seq('-', Bio.Alphabet.generic_dna)
        else:
            out_nucl_recs += [nucl_rec]
    return out_nucl_recs

def codons(nucl_seq):
    """Produces a codon iterator from a nucleotide sequence."""
    for i in range(len(nucl_seq) / 3):
        yield nucl_seq[i * 3:(i * 3) + 3]

def back_align(nucl_rec, aa_rec):
    """Takes a aligned AAs and matches nucleotides to them.
    
    TODO: deal with final stop codons (or lack there of).

    """
    nucl_seq = nucl_rec.seq.ungap(gap = '-')
    aa_seq = aa_rec.seq
    trans_nucl = nucl_seq.translate()
    aa_no_gaps = aa_seq.ungap(gap = '-')
    # remove stop codons from the end of sequences that have them.
    if trans_nucl[-1] == "*":
        nucl_seq = nucl_seq[:-3]
        trans_nucl = nucl_seq.translate()
    if not str(trans_nucl) == str(aa_no_gaps):
        raise ValueError("Nucleotide sequence does not translate to AA-sequence:\nAA (w/out gaps) (%d):\n%s\nTrans-Nucl (%d):\n%s\n" % (len(aa_no_gaps), str(aa_no_gaps), len(trans_nucl), str(trans_nucl)))
    out_seq = Bio.Seq.Seq("", alphabet = Bio.Alphabet.generic_dna)
    codon_gen = codons(nucl_seq)
    for aa in aa_seq:
        if aa is not "-":
            out_seq += codon_gen.next()
        else:
            out_seq += "---"
    return out_seq

def main():
    usage = "usage: %prog [options] [nucl] [aligned-aa]"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option("-n", "--nucl", "--in", dest = "nucl_file",
                      help = "unaligned nucleotides in FASTA format")
    parser.add_option("-a", "--aa", "--alignment", dest = "aa_file",
                      help = "aligned amino acids in FASTA format")
    parser.add_option("-o", "--out", dest = "out_file",
                      help = "filename for aligned nucleotides output")
    (opts, args) = parser.parse_args()
    nucl_file = None
    aa_file = None
    out_file = None
    if opts.nucl_file is not None:
        nucl_file = open(opts.nucl_file)
    elif len(args) >= 1:
        nucl_file = open(args[0])
    else:
        nucl_file = sys.stdin
    if opts.aa_file is not None:
        aa_file = open(opts.aa_file)
    elif len(args) >= 2:
        aa_file = open(args[1])
    else:
        raise TypeError("No aa_file argument or option given")
    if opts.out_file is not None:
        out_file = open(opts.out_file, 'w')
    else:
        out_file = sys.stdout
    nucl_recs = list(Bio.SeqIO.parse(nucl_file, 'fasta', alphabet = Bio.Alphabet.generic_dna))
    aa_recs = list(Bio.SeqIO.parse(aa_file, 'fasta', alphabet = Bio.Alphabet.generic_protein))
    Bio.SeqIO.write(back_align_all(nucl_recs, aa_recs), out_file, 'fasta')

if __name__ == "__main__":
    main()
