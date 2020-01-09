from collections import defaultdict

overlap = defaultdict(int)
t1totals = defaultdict(int)
t2totals = defaultdict(int)
labels = set()
total = 0
for filename in snakemake.input:
    with open(filename) as f:
        for line in f:
            cols = line.split()
            first = cols[3]
            second = cols[7]
            length = int(cols[8])
            overlap[first, second] += length
            overlap[second, first] += length
            total += length
            t1totals[first] += length
            t2totals[first] += length
            t1totals[second] += length
            t2totals[second] += length
            labels.add(first)

with open(snakemake.output.csv, 'w') as outfile:
    outfile.write(' ,' + ','.join(sorted(labels)) + '\n')
    for l1 in sorted(labels):
        line = l1
        for l2 in sorted(labels):
            t1p = t1totals[l1] / float(total)
            t2p = t2totals[l2] / float(total)
            expected = t1p * t2p
            actual = overlap[l1, l2] / float(total)
            line += ',' + str(actual / expected)
        outfile.write(line + '\n')

