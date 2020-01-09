#Prints the locations of every 200 bp segment predicted in the given state
#in at least one of the segmentations given as input

import sys
import collections

state = sys.argv[1]
segmentations = sys.argv[2:]

#Get all positions with state in at least one input segmentation
states = collections.defaultdict(set)
for filename in segmentations:
    with open(filename) as f:
        for line in f:
            cols = line.split()
            if cols[3] == state:
                for num in range(int(cols[1]), int(cols[2]), 200):
                    states[cols[0]].add(num)


for k, v in states.items():
    for num in sorted(v):
        print("{}\t{}\t{}".format(k, num, num+200))

