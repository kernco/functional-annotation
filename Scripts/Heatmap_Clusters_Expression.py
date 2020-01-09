
import math
from collections import defaultdict
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

mapping = {}
exon1 = {}
genelen = {}
startcodon = {}
trnstart = {}
with open(snakemake.input.reference) as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.split()
        if cols[2] == 'transcript':
            gene_id = cols[9].strip(';"')
            trans_id = cols[13].strip(';"')
            mapping[trans_id] = gene_id
            genelen[trans_id] = int(cols[4]) - int(cols[3])
            if cols[6] == '+':
                trnstart[trans_id] = int(cols[3])
            elif cols[6] == '-':
                trnstart[trans_id] = int(cols[4])
        elif cols[2] == 'exon' and cols[17].strip(';"') == '1':
            trans_id = cols[13].strip(';"')
            exon1[trans_id] = int(cols[4]) - int(cols[3])
        elif cols[2] == 'start_codon':
            trans_id = cols[13].strip(';"')
            startcodon[trans_id] = abs(int(cols[3]) - trnstart[trans_id])



refs = {}
with open(snakemake.input.merged) as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.split('\t')
        if cols[2] == 'transcript':
            tags = [x.split() for x in cols[-1].split(';')]
            tags = dict([(x[0], x[1]) for x in tags if x])
            if 'ref_gene_id' in tags:
                refs[tags["ref_gene_id"].strip('"')] = tags["gene_id"].strip('"')

tpms = {}
with open(snakemake.input.expression) as f:
    f.readline()
    for line in f:
        cols = line.split()
        tpms[cols[0]] = float(cols[-1])

clusters = defaultdict(list)
exonlens = defaultdict(list)
genelens = defaultdict(list)
codons = defaultdict(list)
with open(snakemake.input.clusters) as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.split()
        tpm = tpms[refs[mapping[cols[3]]]]
        clusters[cols[-1]].append(math.log(tpm + 1))
        exonlens[cols[-1]].append(exon1[cols[3]])
        genelens[cols[-1]].append(genelen[cols[3]])
        try:
            codons[cols[-1]].append(startcodon[cols[3]])
        except KeyError:
            pass
        except ValueError:
            pass

outfile = open(snakemake.output.txt, 'w')
for k, v in clusters.items():
    outfile.write("{}\t{}\t{}\t{}\t{}\n".format(k, sum(v) / len(v), sum(codons[k]) / len(codons[k]), sum(genelens[k]) / len(genelens[k]), sum(exonlens[k]) / len(exonlens[k])))

def top(lst, cutoff):
    return sorted(lst)[:int(len(lst)*cutoff)]

f, ax = plt.subplots(2,1)
n, bins, patches = ax[0].hist(top(genelens['cluster_1'], 0.9), 50, normed=1, facecolor='blue', alpha=0.5)
#n, bins, patches = ax[0].hist(genelens['cluster_2'], 50, normed=1, facecolor='red', alpha=0.5)
n, bins, patches = ax[0].hist(top(genelens['cluster_3'], 0.9), 50, normed=1, facecolor='green', alpha=0.5)
#n, bins, patches = ax[0].hist(genelens['cluster_4'], 50, normed=1, facecolor='gray', alpha=0.5)
ax[0].set_title('Gene Length')
#plt.savefig('genelens.png')

#f, ax = plt.subplots(2,2)
#n, bins, patches = ax[1].hist(top(exonlens['cluster_1'], 0.9), 50, normed=1, facecolor='blue', alpha=0.5)
#n, bins, patches = ax[1].hist(exonlens['cluster_2'], 50, normed=1, facecolor='red', alpha=0.5)
#n, bins, patches = ax[1].hist(top(exonlens['cluster_3'], 0.9), 50, normed=1, facecolor='green', alpha=0.5)
#n, bins, patches = ax[1].hist(exonlens['cluster_4'], 50, normed=1, facecolor='gray', alpha=0.5)
#ax[1].set_title('Exon 1 Length')
#plt.savefig('exonlens.png')

#f, ax = plt.subplots(2,2)
n, bins, patches = ax[1].hist(top(codons['cluster_1'], 0.75), 50, normed=1, facecolor='blue', alpha=0.5)
#n, bins, patches = ax[2].hist(codons['cluster_2'], 50, normed=1, facecolor='red', alpha=0.5)
n, bins, patches = ax[1].hist(top(codons['cluster_3'], 0.75), 50, normed=1, facecolor='green', alpha=0.5)
#n, bins, patches = ax[2].hist(codons['cluster_4'], 50, normed=1, facecolor='gray', alpha=0.5)
ax[1].set_title('TSS to Start Codon Distance')
plt.savefig('codondist.png')
