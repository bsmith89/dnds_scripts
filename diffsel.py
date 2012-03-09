#!/usr/bin/env python
"""Test differential dN/dS for lineages unique to treatments.

A pipelining of my dN/dS analysis.  Takes a design file, newick
formatted tree, and nucleotide alignment and tests for any difference
between treatments as specified in the design file.

"""
import optparse
import sys
