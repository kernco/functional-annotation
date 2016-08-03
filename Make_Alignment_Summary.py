import sys
import subprocess
import os

print "{: >15s} {: >6s} {: >14s} {: >23s} {: >23s} {: >23s}".format("Tissue", "Rep", "Raw Reads", "Trimmed Reads", "Aligned Reads", "Filtered Reads")
for library in sys.argv[1:]:
    if os.path.isfile("Trimmed_Reads/{}_R2.fq.gz_trimming_report.txt".format(library)):
        reportname = "Trimmed_Reads/{}_R2.fq.gz_trimming_report.txt".format(library)
        readname = "Trimmed_Reads/{}_R2_val_2.fq.gz".format(library)
        divfactor = 2
    elif os.path.isfile("Trimmed_Reads/{}.fq.gz_trimming_report.txt".format(library)):
        reportname = "Trimmed_Reads/{}.fq.gz_trimming_report.txt".format(library)
        readname = "Trimmed_Reads/{}_trimmed.fq.gz".format(library)
        divfactor = 1
    with open(reportname) as f:
        for line in f:
            if line.startswith("Total reads processed"):
                rawreads = int(line.split()[3].replace(',',''))
    zcat = subprocess.Popen(("zcat", readname), stdout=subprocess.PIPE)
    trimreads = int(subprocess.check_output(('wc', '-l'), stdin=zcat.stdout)) / 4
    if library.startswith('RNASeq'):
        alignpairs = int(subprocess.check_output(["samtools", "view", "-c", "-F", "4", "Tophat_Output/{}/accepted_hits.bam".format(library)])) / divfactor
    else:
        alignpairs = int(subprocess.check_output(["samtools", "view", "-c", "-F", "4", "Bwa_Output/{}.aligned.bam".format(library)])) / divfactor
    filterpairs = int(subprocess.check_output(["samtools", "view", "-c", "Aligned_Reads/{}.bam".format(library)])) / divfactor
    print("{: >15s} {: >6s} {: >14,d} {: >14,d} ({: >5.1f}%) {: >14,d} ({: >5.1f}%) {: >14,d} ({: >5.1f}%)".format(library.split('_')[1], library.split('_')[2], rawreads, trimreads, (float(trimreads) / rawreads)*100, alignpairs, (float(alignpairs) / trimreads)*100, filterpairs, (float(filterpairs) / alignpairs)*100))
