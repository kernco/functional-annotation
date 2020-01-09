#!/usr/bin/python env
import collections
import subprocess
import os
import itertools

matrix = collections.defaultdict(dict)
files = list(zip(snakemake.input.peaks, snakemake.input.bedgraphs))
for files1, files2 in itertools.combinations(zip(snakemake.input.peaks, snakemake.input.bedgraphs), 2):
    print(files1, files2)
    peaks1, bedgraph1 = files1
    peaks2, bedgraph2 = files2
    union = subprocess.check_output("cat {} {} | sort -k1,1 -k2,2n | bedtools merge -i stdin | bedtools intersect -a stdin -b {} {} -wb -filenames -sorted | cut -f8 | paste -sd+ | bc".format(peaks1, peaks2, bedgraph1, bedgraph2), shell=True)
    intersection = subprocess.check_output("bedtools intersect -a {} -b {} -sorted | cut -f1,2,3 | bedtools intersect -a stdin -b {} {} -wb -filenames -sorted | cut -f8 | paste -sd+ | bc".format(peaks1, peaks2, bedgraph1, bedgraph2), shell=True)
    jaccard = float(intersection) / float(union)
    if snakemake.params.hide_assay:
        library1 = '_'.join('_'.join(os.path.basename(peaks1).split('.')[:-1]).split('_')[1:])
        library2 = '_'.join('_'.join(os.path.basename(peaks2).split('.')[:-1]).split('_')[1:])
    else:
        library1 = '_'.join('_'.join(os.path.basename(peaks1).split('.')[:-1]).split('_'))
        library2 = '_'.join('_'.join(os.path.basename(peaks2).split('.')[:-1]).split('_'))
    matrix[library1][library2] = jaccard
    matrix[library2][library1] = jaccard
    matrix[library1][library1] = 1
    matrix[library2][library2] = 1

keys = sorted(matrix.keys())
with open(snakemake.output.txt, 'w') as outfile:
    outfile.write("\t" + "\t".join(keys) + '\n')
    for k in keys:
        outfile.write(k)
        for j in keys:
            outfile.write('\t' + matrix[k][j])
        outfile.write('\n')
