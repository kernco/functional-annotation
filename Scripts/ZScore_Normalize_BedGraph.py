import scipy.stats

chromsizes = {}
with open(snakemake.input.chromsizes) as f:
    for line in f:
        cols = line.split()
        chromsizes[cols[0]] = int(cols[1])

segments = []
scores = []
with open(snakemake.input.bedgraph) as f:
    for line in f:
        cols = line.split()
        for start in range(int(cols[1]), int(cols[2]), 100):
            try:
                end = min(start + 100, chromsizes[cols[0]])
                segments.append((cols[0], start, end))
                scores.append(max(0, float(cols[3])))
            except KeyError:
                pass

zscores = scipy.stats.zscore(scores)
with open(snakemake.output.bedgraph, "w") as f:
    for segment, zscore in zip(segments, zscores):
        print('\t'.join([str(x) for x in segment] + [str(zscore)]), file=f)

