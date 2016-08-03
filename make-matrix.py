#!/usr/bin/python env
import sys
import collections 

matrix = collections.defaultdict(dict)
for line in sys.stdin:
	fields = line.strip().split()
	matrix[fields[0]][fields[1]] = float(fields[2])

keys = sorted(matrix.keys())
sys.stdout.write("\t" + "\t".join(keys) + '\n')
for k in keys:
	sys.stdout.write(k)
	for j in keys:
		sys.stdout.write('\t' + str(matrix[k][j]))
        sys.stdout.write('\n')
