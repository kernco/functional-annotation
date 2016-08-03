import sys
import subprocess
import os
import re

samples = sys.argv[1].split()

sys.stdout.write("{: >15s}{: >7s}{: >15s}{: >10s}{: >10s}{: >10s}\n".format("Tissue", "Rep", "Depth", "PBC", "NSC", "RSC"))
for sample in samples:
    output = "{: >15s}{: >7s}".format(*sample.split('_')[1:])
    if os.path.isfile('Aligned_Reads/{}.bam'.format(sample)):
        depth = int(subprocess.check_output(['samtools', 'view', '-c', 'Aligned_Reads/{}.bam'.format(sample)]))
        output += "{: >15,d}".format(depth)
    else:
        output += "{: >15s}".format(" ")
    if os.path.isfile('Bwa_Output/{}.sorted.bam'.format(sample)):
        result = subprocess.check_output(('samtools', 'depth', '-Q 15', 'Bwa_Output/{}.sorted.bam'.format(sample)))
        exact1 = 0
        atleast1 = 0
        for line in result.strip().split('\n'):
            if line:
                cols = line.split()
                if cols[-1] == '1':
                    exact1 += 1
                atleast1 += 1
        output += "{: >10.2f}".format(float(exact1)/atleast1)
    else:
        output += "{: >10s}".format(" ")
    if os.path.isfile('Aligned_Reads/{}.bam'.format(sample)):
        try:
            result = subprocess.check_output(('Rscript', '/home/ckern/phantompeakqualtools/run_spp_nodups.R', '-rf', '-c=Aligned_Reads/{}.bam'.format(sample), '-s=0:2:400', '-savp=Results/{}_Cross_Correlation.pdf'.format(sample), '-p={}'.format(sys.argv[2])))
        except subprocess.CalledProcessError:
            result = ''
        #subprocess.call(('convert', '-geometry 640x640', '-density 200', '-quality 100', 'Results/{}_Cross_Correlation.pdf'.format(sample), 'Results/{}_Cross_Correlation.png'.format(sample)))

        #subprocess.call(('rm', 'Results/{}_Cross_Correlation.pdf'.format(sample)))
        try:
            nsc = float(re.search(r'\(NSC\) ([0-9\.]+)', result).group(1))
            #nsc = float(result.split('\n')[-2].split()[-1])
            output += "{: >10.2f}".format(nsc)
        except (ValueError, IndexError):
            output += "{: >10s}".format('NA')
        try:
            rsc = float(re.search(r'\(RSC\) ([0-9\.]+)', result).group(1))
            #rsc = float(result.split('\n')[-1].split()[-1])
            output += "{: >10.2f}".format(rsc)
        except (ValueError, IndexError):
            output += "{: >10s}".format('NA')
    else:
        output += "{: ^10s}{: ^10s}".format(' ', ' ')
    sys.stdout.write(output + "\n")

