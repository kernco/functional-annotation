import sys
import subprocess
import os

samples = sys.argv[1].split()
assays = sys.argv[2].split()

outfile = open(sys.argv[4], "w")
outfile.write("{: >15s}{: >10s}".format("Tissue", "Rep") + ''.join(["{: ^10s}".format(x) for x in assays]) + "\n")
for sample in samples:
    output = "{: >15s}{: >10s}".format(*sample.split('_'))
    for assay in assays:
        if os.path.isfile('Aligned_Reads/{}_{}.bam'.format(assay, sample)):
            result = subprocess.check_output(('Rscript', '/home/ckern/phantompeakqualtools/run_spp_nodups.R', '-c=Aligned_Reads/{}_{}.bam'.format(assay, sample), '-s=0:2:400', '-p={}'.format(sys.argv[3])))
            score = int(result.split()[-1])
            output += "{: ^ 10d}".format(score)
        else:
            output += "{: ^ 10s}".format(' ')
    outfile.write(output + "\n")

