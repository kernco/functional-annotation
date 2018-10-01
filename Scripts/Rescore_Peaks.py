#Column 5 of the BED format is for a score between 0 and 1000.
#Macs2, however, puts numbers here larger than 1000 which prevents
#formatting the file for use in the UCSC browser. We just set all
#numbers greater than 1000 to 1000.

outfile = open(snakemake.output.scaledbed, 'w')
with open(snakemake.input.narrowpeak) as f:
    for line in f:
        cols = line.split()
        if int(cols[4]) > 1000:
            cols[4] = '1000'
        outfile.write('\t'.join(cols) + '\n')
