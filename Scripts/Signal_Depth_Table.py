from collections import defaultdict
import json

samples = set()
assays = set()
reps = defaultdict(list)
rips = {}
rpgcs = {}

for filename in snakemake.input.stats:
    assay, tissue, rep = filename.split('/')[1].split('_')[:3]
    samples.add(' '.join([tissue, rep]))
    assays.add(assay)
    reps[' '.join([tissue, assay])].append(rep)
    with open(filename) as f:
        stats = json.loads(f.read())
        rips[' '.join([assay, tissue, rep])] = int(stats['Signal Reads'])
        rpgcs[' '.join([assay, tissue, rep])] = float(stats['Signal Reads per KB of Genome'])

txtout = open(snakemake.output.txt, 'w')
csvout = open(snakemake.output.csv, 'w')

headers = ['Tissue', 'Rep'] + sorted(assays)
headerformat = '{: >20s}{: >7s}' + '{: >15s}' * len(assays)
rowformat = '{: >20s}{: >7s}' + '{: >15,d}' * len(assays)
print(headerformat.format(*headers), file=txtout)
print(','.join(headers), file=csvout)
for sample in sorted(samples):
    data = sample.split()
    for assay in sorted(assays):
        try:
            data.append(rips[' '.join([assay, sample])])
        except KeyError:
            data.append(0)
    print(rowformat.format(*data), file=txtout)
    print(','.join([str(x) for x in data]), file=csvout)

print('\n' + headerformat.format(*headers), file=txtout)
print('\n' + ','.join(headers), file=csvout)
rowformat = '{: >20s}{: >7s}' + '{: >15.1f}' * len(assays)
for sample in sorted(samples):
    data = sample.split()
    for assay in sorted(assays):
        try:
            data.append(rpgcs[' '.join([assay, sample])])
        except KeyError:
            data.append(0)
    print(rowformat.format(*data), file=txtout)
    print(','.join([str(x) for x in data]), file=csvout)

