import sys

badids = set()
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('#'):
            print line.strip()
            continue
        if 'C_gene_segment' in line:
            badids.add(line.split('\t')[-1].split(';')[0].split('=')[1])
        elif line.split('\t')[-1].split(';')[1].split('=')[1] in badids:
            continue
        else:
            print line.strip()

