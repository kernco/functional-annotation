from collections import defaultdict
import json

samples = set()
assays = set()
reps = defaultdict(list)
rips = {}
rpgcs = {}

for filename in set(snakemake.input.stats):
    assay, tissue = filename.split('/')[1].split('_')[:2]
    samples.add(tissue)
    assays.add(assay)
    #reps[' '.join([tissue, assay])].append(rep)
    with open(filename) as f:
        f.readline() #Discard headers
        rip = 0
        total = 0
        for line in f:
            cols = line.split()
            rip += int(cols[-2])
            total += int(cols[-1])
        rips[' '.join([assay, tissue])] = rip
        rpgcs[' '.join([assay, tissue])] = (rip * 1000) / snakemake.config['genomesize']

txtout = open(snakemake.output.txt, 'w')
csvout = open(snakemake.output.csv, 'w')

headers = ['Tissue'] + sorted(assays)
headerformat = '{: >20s}' + '{: >15s}' * len(assays)
rowformat = '{: >20s}' + '{: >15,d}' * len(assays)
print(headerformat.format(*headers), file=txtout)
print(','.join(headers), file=csvout)
for sample in sorted(samples):
    data = [sample]
    for assay in sorted(assays):
        try:
            data.append(rips[' '.join([assay, sample])])
        except KeyError:
            data.append(0)
    print(rowformat.format(*data), file=txtout)
    print(','.join([str(x) for x in data]), file=csvout)

print('\n' + headerformat.format(*headers), file=txtout)
print('\n' + ','.join(headers), file=csvout)
rowformat = '{: >20s}' + '{: >15.1f}' * len(assays)
for sample in sorted(samples):
    data = [sample]
    for assay in sorted(assays):
        try:
            data.append(rpgcs[' '.join([assay, sample])])
        except KeyError:
            data.append(0)
    print(rowformat.format(*data), file=txtout)
    print(','.join([str(x) for x in data]), file=csvout)

