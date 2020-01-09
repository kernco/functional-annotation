import subprocess

def count_lines(filename):
    output = subprocess.check_output(['wc', '-l', filename])
    return int(output.split()[0])

def count_merged_regions(filename):
    output = subprocess.check_output(['bedtools', 'merge', '-i', filename])
    return len(output.decode('utf-8').split('\n'))

def sum_region_sizes(filename):
    size = 0
    with open(filename) as f:
        for line in f:
            cols = line.split()
            size += int(cols[2]) - int(cols[1])
    return size

gsize = float(snakemake.config['genomesize'])
csvout = open(snakemake.output[0], 'w')
for tissue in sorted(snakemake.params.tissues):
    cols = [tissue]
    for rep in sorted(snakemake.params.reps):
        repid = '{}_{}_{}'.format(snakemake.wildcards.assay, tissue, rep)
        peakfile = 'Enriched_Regions/{}.bed'.format(repid)
        regions = count_lines(peakfile)
        #cols.append(str(regions))
        cols.append(str(count_merged_regions(peakfile)))
        cols.append(str((regions*200) / gsize))
    peakfile = 'Enriched_Regions/{}_{}_Combined.bed'.format(snakemake.wildcards.assay, tissue)
    regions = count_lines(peakfile)
    cols.append(str(regions))
    #cols.append(str(count_merged_regions(peakfile)))
    cols.append(str(sum_region_sizes(peakfile) / gsize))
    csvout.write(','.join(cols) + '\n')
