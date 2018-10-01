#!/usr/bin/python env
import collections
import subprocess
import os

matrix = collections.defaultdict(dict)
for peaks1 in snakemake.input.peaks:
    for peaks2 in snakemake.input.peaks:
        result = subprocess.getoutput("bedtools jaccard -a {} -b {}".format(peaks1, peaks2))
        try:
            jaccard = result.split('\n')[1].split()[2]
        except IndexError:
            jaccard = '0'
        library1 = '_'.join('_'.join(os.path.basename(peaks1).split('.')[:-1]).split('_')[1:])
        library2 = '_'.join('_'.join(os.path.basename(peaks2).split('.')[:-1]).split('_')[1:])
        matrix[library1][library2] = jaccard

keys = sorted(matrix.keys())
with open(snakemake.output.txt, 'w') as outfile:
    outfile.write("\t" + "\t".join(keys) + '\n')
    for k in keys:
            outfile.write(k)
            for j in keys:
                    outfile.write('\t' + matrix[k][j])
            outfile.write('\n')
