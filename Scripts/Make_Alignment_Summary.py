import json

txtout = open(snakemake.output.txt, 'w')
csvout = open(snakemake.output.csv, 'w')
txtout.write("{: >15s} {: >6s} {: >14s} {: >23s} {: >23s} {: >23s} {: >23s} {: >23s}\n".format("Tissue", "Rep", "Raw Reads", "Trimmed", "Unaligned", "Filtered", "Duplicates", "Remaining"))
headers = ["Tissue", "Rep", "Raw Reads", "Trimmed", "%", "Unaligned", "%", "Filtered", "%", "Duplicates", "%", "Remaining", "%"]
csvout.write(','.join(headers) + '\n')
for library in snakemake.params.libraries:
    assay, tissue, replicate = library.split('_')
    with open('Metrics/{}_Alignment_Stats.json'.format(library)) as f:
        stats = json.loads(f.read())
    rawreads = int(stats['Raw Reads'])
    trimreads = int(stats['Trimmed Reads'])
    alignpairs = int(stats['Aligned Reads'])
    filterpairs = int(stats['Filtered Reads'])
    try:
        deduped = int(stats['Deduplicated Reads'])
    except KeyError:
        deduped = filterpairs
    data = [tissue, replicate, rawreads, rawreads - trimreads, (float(rawreads - trimreads) / rawreads)*100, trimreads - alignpairs, (float(trimreads - alignpairs) / trimreads)*100, alignpairs - filterpairs, (float(alignpairs - filterpairs) / alignpairs)*100, filterpairs - deduped, (float(filterpairs - deduped) / filterpairs)*100, deduped, (float(deduped) / rawreads)*100]
    txtout.write("{: >15s} {: >6s} {: >14,d} {: >14,d} ({: >5.1f}%) {: >14,d} ({: >5.1f}%) {: >14,d} ({: >5.1f}%) {: >14,d} ({: >5.1f}%) {: >14,d} ({: >5.1f}%)\n".format(*data))
    csvout.write(','.join([str(x) for x in data]) + '\n')
