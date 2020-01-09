import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn import metrics
import sys
#import math

def load_data(filename):
    data = []
    with open(filename) as f:
        for line in f:
            cols = line.strip().split('\t')
            data.append((float(cols[0]), 1 if cols[1] == "True" else 0))
    return data

hcp_k4me3 = load_data(sys.argv[1])
hcp_k27ac = load_data(sys.argv[2])
lcp_k4me3 = load_data(sys.argv[3])
lcp_k27ac = load_data(sys.argv[4])

def plot_curve(data, name):
    fpr, tpr, thresholds = metrics.roc_curve([x[1] for x in data], [x[0] for x in data], pos_label=1)
    auc = metrics.auc(fpr, tpr)
    #predscore += len(gene_lists[state]) * abs(auc - 0.5) * 2
    #predtotal += len(gene_lists[state])
    plt.plot(fpr, tpr, label='{} (area = {:0.2f})'.format(name, auc))

#ROC
fig = plt.figure(figsize=(8,8), dpi=100)
#predscore = 0
#predtotal = 0
plot_curve(hcp_k4me3, "HCPs w/H3K4me3")
plot_curve(hcp_k27ac, "HCPs w/H3K27ac")
plot_curve(lcp_k4me3, "LCPs w/H3K4me3")
plot_curve(lcp_k27ac, "LCPs w/H3K27ac")
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(loc="best", prop={'size': 10})
#plt.title("{}".format(predscore / predtotal))
fig.savefig('test.png', bbox_inches='tight')
plt.close(fig)

