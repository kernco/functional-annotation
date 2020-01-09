from collections import defaultdict

outfile = open(snakemake.output.tss, 'w')
tss = defaultdict(list)
with open(snakemake.input.gtf) as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.strip().split('\t')
        try:
            if cols[2] == 'transcript':
                gene_id = cols[-1].split(';')[0].split()[1].strip('"')
                transcript_id = cols[-1].split(';')[1].split()[1].strip('"')
                if cols[6] == '+':
                    tss[gene_id + " " + cols[6] + ' ' + cols[0]].append(int(cols[3]))
                else:
                    tss[gene_id + " " + cols[6] + ' ' + cols[0]].append(int(cols[4]))
        except IndexError:
            pass

with open(snakemake.output.tss, 'w') as f:
    for k, v in tss.items():
        gene_id, strand, chrom = k.split()
        if strand == '-':
            coord = max(v)
        else:
            coord = min(v)
        outfile.write('\t'.join([chrom, str(coord), str(coord + 1), gene_id, '0', strand]) + '\n')

