import sys
import re

with open(sys.argv[1]) as f:
    for line in f:
        parts = line.split('\t')
        gene_id = re.search(r'gene_id "(XLOC_[0-9]+)"', line).group(1)
        trans_id = re.search(r'transcript_id "(TCONS_[0-9]+)"', line).group(1)
        exon_num = re.search(r'exon_number "([0-9]+)"', line).group(1)
        name = "{}.{}.exon{}".format(gene_id, trans_id, exon_num)
        if parts[3] != parts[4]:
            print "\t".join([parts[0], parts[3], parts[4], name, '0', parts[6]])
