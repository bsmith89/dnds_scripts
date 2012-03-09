#!/usr/bin/env python
"""Randomize sample-treatment mapping.

Takes space delimitted design file to stdin and prints a tab delimited,
randomized design file to stdout.

Meant for making sure that the null hypothesis holds for de-associated
data.

"""

import sys
import random

samples = []
groups = []
for mapping in sys.stdin:
    sample, group = mapping.split()
    samples += [sample]
    groups += [group]
# now reassign all of the groupings
rand_groups = random.sample(groups, len(groups))
for sample, group in zip(samples, rand_groups):
    sys.stdout.write("%s\t%s\n" % (sample, group))
