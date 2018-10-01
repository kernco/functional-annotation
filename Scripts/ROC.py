import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn import metrics
import sys
#import math

def load_tpms(filename):
    tpms = {}
    with open(filename) as f:
        f.readline()
        for line in f:
            cols = line.strip().split()
            tpms[cols[0]] = float(cols[-1])
    return tpms

reps = []
for filename in sys.argv[2:]:
    reps.append(load_tpms(filename))

gene_lists = defaultdict(list)
with open(sys.argv[1]) as f:
    for line in f:
        cols = line.split('\t')
        gene_id = cols[3]
        state = int(cols[9])
        gene_lists[state].append(gene_id)

#ROC
fig = plt.figure(figsize=(8,8), dpi=100)
predscore = 0
predtotal = 0
for state in reversed(sorted(gene_lists)):
    data = []
    for tpm in reps:
        for gene_id in tpm:
            if tpm[gene_id] > 0:
                data.append((tpm[gene_id], 1 if gene_id in gene_lists[state] else 0))
    scores, preds = zip(*sorted(data))
    fpr, tpr, thresholds = metrics.roc_curve(preds, scores, pos_label=1)
    auc = metrics.auc(fpr, tpr)
    predscore += len(gene_lists[state]) * abs(auc - 0.5) * 2
    predtotal += len(gene_lists[state])
    plt.plot(fpr, tpr, label='State {} (area = {:0.2f})'.format(state, auc))
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(loc="best", prop={'size': 10})
plt.title("{}".format(predscore / predtotal))
fig.savefig('test.png', bbox_inches='tight')
plt.close(fig)

