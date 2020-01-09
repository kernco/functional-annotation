import sys
import os
from collections import defaultdict

def load_count_table(filename):
    groups = defaultdict(set)
    higher = defaultdict(set)
    with open(filename) as f:
        headers = [x.split('_')[0] for x in f.readline().split()]
        for line in f:
            cols = line.split()
            cpms = list(reversed(sorted([(float(x),tissue) for x, tissue in zip(cols[1:], headers)])))
            if cpms[0][0] >= 1 and cpms[1][0]*1.5 < cpms[0][0]:
                higher[cpms[0][1]].add(cols[0])
            #Tissue enriched
            if cpms[0][0] >= 1 and cpms[1][0]*4 < cpms[0][0]:
                groups[cpms[0][1]].add(cols[0])
            #Group enriched
            else:
                for i in range(2,len(cpms)):
                    for cpm in cpms[:i]:
                        if cpm[0] >= 1 and cpms[i][0]*4 > cpm[0]:
                            break
                    else:
                        groups['_'.join(sorted([x[1] for x in cpms[:i]]))].add(cols[0])
                        break
    return groups, higher

groups, higher = load_count_table(snakemake.input.counts)

os.system("mkdir -p {}".format(snakemake.output.outdir))
for group, genes in groups.items():
    with open(snakemake.output.outdir + '/' + group + '_Enriched_Genes.txt', 'w') as groupout:
        print('\n'.join(genes), file=groupout)

for group, genes in higher.items():
    with open(snakemake.output.outdir + '/' + group + '_Higher_Genes.txt', 'w') as groupout:
        print('\n'.join(genes), file=groupout)
