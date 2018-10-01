
refs = {}
with open(snakemake.input.gtf) as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.split('\t')
        if cols[2] == 'transcript':
            tags = [x.split() for x in cols[-1].split(';')]
            tags = dict([(x[0], x[1]) for x in tags if x])
            if 'ref_gene_id' in tags:
                refs[tags["gene_id"].strip('"')] = tags["ref_gene_id"].strip('"')
            else:
                refs[tags["gene_id"].strip('"')] = '-'


tpms = {}
genes = set()
for filename in snakemake.input.tabs:
    with open(filename) as f:
        f.readline()
        for line in f:
            cols = line.split()
            tpms[cols[0], filename] = cols[-1]
            genes.add(cols[0])

with open(snakemake.output.tsv, 'w') as outfile:
    outfile.write("GeneID\tRefID\t")
    for filename in snakemake.input.tabs:
        header = "_".join(filename.split("_")[1:-1])
        outfile.write(header + "\t")
    outfile.write("\n")
    for gene in sorted(genes):
        outfile.write("{}\t{}\t".format(gene, refs[gene]))
        for filename in snakemake.input.tabs:
            outfile.write(tpms[gene, filename] + "\t")
        outfile.write("\n")
outfile.close()

