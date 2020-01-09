import seaborn
import math
import matplotlib.pyplot as plt
from collections import defaultdict

colors = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#991eb4", "#42d4f4", "#f032e6"]

def make_plot(filenames, outfile, xaxis_label):
    for filename in filenames:
        data = defaultdict(list)
        with open(filename) as f:
            headers = f.readline().split()[2:]
            #headers = f.readline().split()
            for line in f:
                cols = line.split()
                #for tissue, tpm in zip(headers, cols[1:]):
                for tissue, tpm in zip(headers, cols[2:]):
                    if float(tpm) > 0:
                        #data[tissue].append(math.log2(float(tpm)+1))
                        data[tissue].append(float(tpm))
        for tissue, color in zip(headers, colors):
            if filename == filenames[0]:
                label = tissue.split('_')[0]
            else:
                label = None
            kplt = seaborn.kdeplot(data[tissue], label=label, color=color)
            kplt.set(xlim=(0,None))

    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    plt.legend()
    #plt.xlabel("TMM-normalized CPM (log2)")
    plt.xlabel(xaxis_label)
    plt.ylabel("Density")
    plt.savefig(outfile, bbox_inches="tight", dpi=400)

#if __name__ == "__main__":
    #import sys
    #make_plot(sys.argv[3:], sys.argv[1], sys.argv[2])
#else:
make_plot(snakemake.input.tpms, snakemake.output.png, "TMM-normalized CPM (log2)")

