import gzip
import sys
import json
from collections import defaultdict
from statistics import mean, StatisticsError
import numpy

#counts = {}
#with open(snakemake.input.chromlens) as f:
    #for line in f:
        #cols = line.split()
        #length = int(cols[1])
        #if length % 200 == 0:
            #length -= 1 #Avoid extra bin when chromosome ends on bin boundary
        #for pos in range(0, length // 200 + 1):
            #counts[(cols[0], pos)] = 0

#for tagalign, stats in zip(snakemake.input.reads, snakemake.input.stats):
    #with open(stats) as f:
        #data = json.loads(f.read())
        #offset = int(data["Est. Fragment Length"]) // 2

    #with gzip.open(tagalign) as f:
        #for read in f:
            #cols = read.decode('utf-8').split()
            #pos = int(cols[1]) + offset
            #try:
                #counts[(cols[0], pos // 200)] += 1
            #except KeyError:
                #pass #Read outside chromosome end

#with open(snakemake.output.counts, 'w') as f:
    #for pos, count in counts.items():
        #print('\t'.join([pos[0], str(pos[1] * 200), str(pos[1] * 200 + 200), str(count)]), file=f)


#With smoothing
counts = defaultdict(list)
chromlens = {}
with open(snakemake.input.chromlens) as f:
    for line in f:
        cols = line.split()
        length = int(cols[1])
        chromlens[cols[0]] = length
        if length % 50 == 0:
            length -= 1 #Avoid extra bin when chromosome ends on bin boundary
        for pos in range(0, length // 50):
            counts[cols[0]].append(0)

for tagalign, stats in zip(snakemake.input.reads, snakemake.input.stats):
    with open(stats) as f:
        data = json.loads(f.read())
        offset = int(data["Est. Fragment Length"]) // 2

    with gzip.open(tagalign) as f:
        for read in f:
            cols = read.decode('utf-8').split()
            pos = int(cols[1]) + offset
            try:
                if cols[0] in counts:
                    counts[cols[0]][pos // 50] += 1
            except IndexError:
                pass #Read outside chromosome end
            #try:
                #counts[cols[0]][pos // 50 + 1] += 0.25
            #except IndexError:
                #pass
            #try:
                #counts[cols[0]][pos // 50 - 1] += 0.25
            #except IndexError:
                #pass

with open(snakemake.output.counts, 'w') as f:
    for chrom, nums in counts.items():
        smoothednums = numpy.convolve(nums, [1/4, 1/2, 1/4], mode='same')
        #try:
            #smoothednums = [mean(nums[:2])] + [mean(nums[i-1:i+2]) for i in range(1,len(nums))]
        #except StatisticsError:
            #print(chrom)
            #sys.exit(1)
        smoothednums = numpy.append(smoothednums, [0] * (4 - len(smoothednums) % 4)) #Pad the end of the list with 0s to a multiple of 4
        bincounts = numpy.sum(numpy.reshape(smoothednums, (-1,4)), axis=1)
        #if chromlens[chrom] % 200 == 0:
            #bincounts = bincounts[:-1]
        for i, count in enumerate(bincounts):
            print('\t'.join([chrom, str(i * 200), str(i * 200 + 200), str(int(round(count)))]), file=f)




