outfile = open(snakemake.output.hets, 'w')
rsid = 1
print('ID\tCHROM\tPOS\tREF\tALT', file=outfile)
with open(snakemake.input.vcf) as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.split()
        print('{}\t{}\t{}\t{}\t{}'.format(str(rsid), cols[0], cols[1], cols[4], cols[5]), file=outfile)
        rsid += 1

