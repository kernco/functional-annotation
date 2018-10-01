
lens = {}
with open(snakemake.input.chromsizes) as f:
    for line in f:
        cols = line.split()
        lens[cols[0]] = int(cols[1])

outfile = open(snakemake.output.trimmed, "w")
with open(snakemake.input.bed) as f:
    for line in f:
        cols = line.split()
        if int(cols[2]) > lens[cols[0]]:
            cols[2] = str(lens[cols[0]])
        outfile.write('\t'.join(cols) + "\n")
