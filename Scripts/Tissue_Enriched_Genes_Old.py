import sys
import os
from collections import defaultdict

def load_count_table(filename):
    groups = defaultdict(set)
    with open(filename) as f:
        headers = [x.split('_')[0] for x in f.readline().split()]
        for line in f:
            cols = line.split()
            cpms = list(reversed(sorted([(float(x),tissue) for x, tissue in zip(cols[1:], headers)])))
            #Tissue enriched
            if cpms[0][0] >= 1 and cpms[1][0]*4 < cpms[0][0]:
                groups[cpms[0][1]].add(cols[0])
            #Group enriched
            else:
                for i in range(2,len(cpms)):
                    for cpm in cpms[:i]:
                        if cpms[i][0]*4 > cpm[0]:
                            break
                    else:
                        groups['_'.join(sorted([x[1] for x in cpms[:i]]))].add(cols[0])
                        break
    return groups

rep1 = load_count_table(snakemake.input.reps[0])
rep2 = load_count_table(snakemake.input.reps[1])

os.system("mkdir -p {}".format(snakemake.output.outdir))
lines = []
for k in set(rep1) | set(rep2):
    intersect = rep1[k] & rep2[k]
    jaccard = len(intersect) / ((len(rep1[k]) + len(rep2[k])) - len(intersect))
    lines.append((len(intersect), jaccard, len(rep1[k]), len(rep2[k]), k))
    if len(intersect) > 0:
        groupout = open(snakemake.output.outdir + '/' + k + '_Enriched_Genes.txt', 'w')
        for gene in intersect:
            print(gene, file=groupout)
        groupout.close()

with open(snakemake.output.summary, 'w') as outfile:
    for line in reversed(sorted(lines)):
        print('{}\t{:.2f}\t{}\t{}\t{}'.format(*line), file=outfile)
