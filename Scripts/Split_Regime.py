
lines = []
with open(snakemake.input.infile) as f:
    f.readline()
    for line in f:
        cols = line.strip().split('\t')
        cols[1] = int(float(cols[1]))
        cols[2] = int(float(cols[2])) #Gets rid of scientific notation, keep int for later sorting
        if cols[3].startswith('regimeR_{}'.format(snakemake.wildcards.regime)):
            lines.append(tuple(cols))

with open(snakemake.output.outfile, 'w') as f:
    for line in sorted(lines):
        print('\t'.join([str(x) for x in line[:3]]), file=f)

