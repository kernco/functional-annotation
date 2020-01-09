import seaborn
import math
import matplotlib.pyplot as plt
from collections import defaultdict

def make_plot(filenames, png):
    for filename in filenames:
        data = []
        with open(filename) as f:
            headers = f.readline()
            for line in f:
                cols = line.split()
                tpm = float(cols[-2])
                tsi = float(cols[0])
                if tpm > 0:
                    data.append(tsi)

    kplt = seaborn.kdeplot(data)
    kplt.set(xlim=(0,1))
    plt.xlabel("Tissue-specificity (Tau)")
    plt.ylabel("Density")
    plt.savefig(png, bbox_inches="tight")

#if __name__ == "__main__":
    #import sys
    #make_plot(sys.argv[1:-1], sys.argv[-1])
#else:
make_plot(snakemake.input.tsis, snakemake.output.png)
