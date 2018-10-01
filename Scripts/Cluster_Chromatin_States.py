import sys
import collections
import numpy
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

#arg1 = state locations
#arg2 = state column (0 indexed)
#arg3+ = tissue names

segments = collections.defaultdict(list)
with open(snakemake.input.states) as f:
    for line in f:
        cols = line.split()
        segments[cols[0]].append(cols[1])

data = []
locs = []
states = snakemake.wildcards.state.split(',')
for k, v in segments.items():
    scores = collections.defaultdict(list)
    for tissue in sorted(snakemake.config['tissues']):
        try:
            with open('ChromHMM/Model_Joint_Reads_{states}/POSTERIOR/{tissue}_{states}_{chrom}_posterior.txt'.format(tissue=tissue, states=snakemake.wildcards.states, chrom=k)) as f:
                f.readline()
                f.readline()
                for line in f:
                    score = 0
                    cols = line.split()
                    for state in states:
                        state = int(state[1:]) - 1 #Remove leading letter
                        score += float(cols[state])
                    scores[tissue].append(score)
        except IOError:
            break
    else:
        for segment in v:
            data.append([scores[tissue][int(segment) // 200] for tissue in sorted(snakemake.config['tissues'])])
            locs.append([k, segment])

array = numpy.array(data)
kmeans = KMeans(n_clusters=16, random_state=0).fit(array)
clusters = collections.defaultdict(list)
outfile = open(snakemake.output.txt, 'w')
for label, scores, loc in zip(kmeans.labels_, data, locs):
    clusters[label].append(scores)
    outfile.write('\t'.join([loc[0], loc[1], str(int(loc[1]) + 1)] + [str(label)]) + '\n')

f, axarr = plt.subplots(figsize=(4,8), nrows=len(clusters), gridspec_kw=dict(height_ratios=[len(c) for c in clusters.values()]))
plt.yticks(visible=False)
for cluster, v in clusters.items():
    axarr[cluster].pcolor(v, cmap=plt.cm.Reds)
    axarr[cluster].set_xticks([])
    axarr[cluster].set_yticks([])
    #axarr[cluster].set_ylabel(str(cluster), rotation=0, size='large', va='bottom')
    axarr[cluster].annotate(chr(65 + cluster), xy=(0, 0.5), xytext=(-5, 0), xycoords='axes fraction', textcoords='offset points', va='center', ha='right')
axarr[-1].set_xticks([x + 0.5 for x in range(len(snakemake.config['tissues']))])
axarr[-1].set_xticklabels(sorted(snakemake.config['tissues']))
#plt.xticks([x + 0.5 for x in range(len(sys.argv[3:]))], sys.argv[3:])
plt.savefig(snakemake.output.png, bbox_inches='tight')
