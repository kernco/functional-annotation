import sys
from collections import defaultdict

labels = {}
for label in sys.argv[3:]:
    parts = label.split(':')
    labels[parts[1]] = parts[0]

colors = {}
with open(sys.argv[2]) as f:
    header = f.readline().strip().split('\t')[1:]
    for line in f:
        parts = line.split()
        cols = {}
        for i,n in enumerate(parts[1:]):
            cols[header[i]] = int(float(n)*255)
        try:
            colors[parts[0]] = "{},{},{}".format(cols[labels['red']], cols[labels['green']], cols[labels['blue']])
        except KeyError:
            print header
            print cols
            print labels
            sys.exit(-1)

lines = defaultdict(list)
with open(sys.argv[1]) as f:
    header = f.readline()
    for line in f:
        parts = line.split('_')
        tissue = parts[0]
        newline = '_'.join(parts[1:])
        lines[tissue].append(newline)

for k, v in lines.iteritems():
    with open('{}.bed'.format(k), 'w') as f:
        f.write('track name="{0}" description="{0}" visibility=1 itemRgb="On"\n'.format(k))
        for line in v:
            parts = line.split('\t')
            parts[-1] = colors[parts[3]]
            f.write('\t'.join(parts) + '\n')

