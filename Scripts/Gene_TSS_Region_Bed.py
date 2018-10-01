data = []
with open(snakemake.input.tab) as f:
    f.readline()
    for line in f:
        cols = line.split()
        if cols[3] == '+':
            tss = cols[4]
        elif cols[3] == '-':
            tss = cols[5]
        else:
            continue
        #if '_' in cols[2]:
            #continue
        data.append((float(cols[-1]), cols[2], tss))

chromsizes = {}
with open(snakemake.input.chromsizes) as f:
    for line in f:
        cols = line.split()
        chromsizes[cols[0]] = int(cols[1])

with open(snakemake.output.bed, 'w') as outfile:
    for tpm, chrom, tss in reversed(sorted(data)):
        if int(tss) - 3000 <= 0 or int(tss) + 3000 >= chromsizes[chrom]:
            continue
        outfile.write("{}\t{}\t{}\n".format(chrom, tss, int(tss) + 1))

