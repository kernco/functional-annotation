tss = []
mapping = {}
with open(snakemake.input.annotation) as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.split('\t')
        tags = {k:v for k, v in [x.strip('"').split(' "')[:2] for x in cols[-1].strip(';').split('; ')]}
        if cols[2] == 'transcript' and 'gene_name' in tags:
            if cols[6] == '+':
                tss.append([cols[0], cols[3], tags['transcript_id']])
            elif cols[6] == '-':
                tss.append([cols[0], cols[4], tags['transcript_id']])
            mapping[tags['transcript_id']] = tags['gene_name']

outfile = open(snakemake.output.gtf, 'w')
for chrom, loc, trn in tss:
    try:
        outfile.write("{}\t{}\t{}\t{}\n".format(chrom, loc, int(loc) + 1, mapping[trn]))
    except KeyError:
        pass



