import sys

with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('#'):
            continue
        elif not line.startswith('NC_'):
            continue
        cols = line.split('\t')
        if cols[2] == 'gene':
            gene_id = cols[-1].split(';')[0].split('=')[1].strip('"')
            sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(cols[0], cols[3], cols[4], 'NA', 0, cols[6]))

