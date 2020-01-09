
seqs = {}
with open(snakemake.input.genome) as f:
    for line in f:
        if line.startswith(">"):
            name = line.split()[0][1:]
            seqs[name] = 0
        else:
            seqs[name] += len(line.strip())

with open(snakemake.output.chroms, 'w') as f:
    for seq in seqs:
        if seq != 'chrM':
            f.write("{}\t{}\n".format(seq, seqs[seq]))
