import matplotlib.pyplot as plt
from collections import defaultdict
#import math

def load_tpms(filename, tpms):
    #tpms = defaultdict(list)
    with open(filename) as f:
        f.readline()
        for line in f:
            cols = line.strip().split()
            tpms[cols[0]].append(float(cols[-1]))
    #return tpms

#reps = []
tpms = defaultdict(list)
for filename in snakemake.input.expression:
    load_tpms(filename, tpms)

gene_lists = defaultdict(list)
with open(snakemake.input.overlap) as f:
    for line in f:
        cols = line.split('\t')
        gene_id = cols[3]
        state = int(cols[9])
        gene_lists[state].append(gene_id)

fig = plt.figure(figsize=(8,8), dpi=100)
data = [[sum(tpms[gene]) / len(tpms[gene]) for gene in gene_lists[state]] for state in reversed(sorted(gene_lists))]
bp = plt.boxplot(data, sym='', vert=False, widths=0.6)
labels = ['State {}: {}'.format(x, len(gene_lists[x])) for x in reversed(sorted(gene_lists))]
ax = plt.axes()
ax.set_yticklabels(labels)
plt.xlabel('Expression (TPM)')
fig.savefig(snakemake.output.png, bbox_inches='tight')
plt.close(fig)

