import subprocess
import os

txtout = open(snakemake.output.txt, 'w')
csvout = open(snakemake.output.csv, 'w')
txtout.write("{: >15s} {: >6s} {: >14s} {: >23s} {: >23s} {: >23s} {: >23s} {: >23s}\n".format("Tissue", "Rep", "Raw Reads", "Trimmed", "Unaligned", "Filtered", "Duplicates", "Remaining"))
headers = ["Tissue", "Rep", "Raw Reads", "Trimmed", "%", "Unaligned", "%", "Filtered", "%", "Duplicates", "%", "Remaining", "%"]
csvout.write(','.join(headers) + '\n')
for library in snakemake.params.libraries:
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
    trimreads = int(int(subprocess.check_output(('wc', '-l'), stdin=zcat.stdout)) / 4)
    if library.startswith('RNASeq'):
        uniquely = subprocess.check_output(["grep", "Uniquely mapped reads number", "STAR_Output/{}_Log.final.out".format(library)])
        multimap = subprocess.check_output(["grep", "Number of reads mapped to multiple loci", "STAR_Output/{}_Log.final.out".format(library)])
        alignpairs = int(uniquely.split(b'|')[1]) + int(multimap.split(b'|')[1])
        deduped = int(int(subprocess.check_output(["samtools", "view", "-c", "Aligned_Reads/{}.bam".format(library)])))
    else:
        alignpairs = int(int(subprocess.check_output(["samtools", "view", "-c", "-F", "4", "Bwa_Output/{}.aligned.bam".format(library)])) / divfactor)
        filterpairs = int(int(subprocess.check_output(["samtools", "view", "-c", "-F", "4", "Bwa_Output/{}.filtered.bam".format(library)])) / divfactor)
        deduped = int(int(subprocess.check_output(["samtools", "view", "-c", "Aligned_Reads/{}.bam".format(library)])) / divfactor)
    tissue = library.split('_')[1]
    try:
        replicate = library.split('_')[2]
    except IndexError:
        replicate = ' '
    if library.startswith('RNASeq'):
        data = [tissue, replicate, rawreads, rawreads - trimreads, (float(rawreads - trimreads) / rawreads)*100, rawreads - alignpairs, (float(rawreads - alignpairs) / rawreads)*100, rawreads - deduped, (float(rawreads - deduped) / rawreads)*100, 0, 0, deduped, (float(deduped) / rawreads)*100]
    else:
        data = [tissue, replicate, rawreads, rawreads - trimreads, (float(rawreads - trimreads) / rawreads)*100, rawreads - alignpairs, (float(rawreads - alignpairs) / rawreads)*100, rawreads - filterpairs, (float(rawreads - filterpairs) / rawreads)*100, rawreads - deduped, (float(rawreads - deduped) / rawreads)*100, deduped, (float(deduped) / rawreads)*100]
    txtout.write("{: >15s} {: >6s} {: >14,d} {: >14,d} ({: >5.1f}%) {: >14,d} ({: >5.1f}%) {: >14,d} ({: >5.1f}%) {: >14,d} ({: >5.1f}%) {: >14,d} ({: >5.1f}%)\n".format(*data))
    csvout.write(','.join([str(x) for x in data]) + '\n')
