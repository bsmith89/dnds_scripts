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
