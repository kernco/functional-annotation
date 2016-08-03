import sys

expressed = set()
notexpressed = set()
with open(sys.argv[1]) as f:
    f.readline()
    for line in f:
        cols = line.split('\t')
        if float(cols[9]) >= float(sys.argv[3]):
            expressed.add(cols[0])
        elif float(cols[9]) <= float(sys.argv[4]):
            notexpressed.add(cols[0])

expfile = open("{}/Expressed_Genes.bed".format(sys.argv[5]), "w")
notfile = open("{}/Repressed_Genes.bed".format(sys.argv[5]), "w")
with open(sys.argv[2]) as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.split('\t')
        if cols[2] == 'gene':
            gene_id = cols[-1].split(';')[0].split('=')[1].strip('"')
            if gene_id in expressed:
                expfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(cols[0], cols[3], cols[4], 'NA', 0, cols[6]))
            elif gene_id in notexpressed:
                notfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(cols[0], cols[3], cols[4], 'NA', 0, cols[6]))

expfile.close()
notfile.close()
