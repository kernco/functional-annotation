import seaborn
import matplotlib.pyplot as plt
from collections import defaultdict
import gzip

colors = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#991eb4", "#42d4f4", "#f032e6"]
def make_graphs(inputfile, depthpng, methpng):
    depths = defaultdict(list)
    methpct = defaultdict(list)
    with gzip.open(inputfile) as f:
        for line in f:
            cols = line.split()
            #depths[cols[3]].append(int(cols[-1]))
            if int(cols[-1]) >= 10:
                methpct[cols[3]].append(float(cols[-3]))
    #for group, color in zip(depths.keys(), colors):
        #kplt = seaborn.kdeplot(depths[group], label=group, color=color)
        #kplt.set(xlim=(0,None))
    #plt.legend()
    #plt.xlabel("Read depth")
    #plt.ylabel("Density")
    #plt.savefig(depthpng, bbox_inches="tight", dpi=400)
    #plt.close()
    for group, color in zip(methpct.keys(), colors):
        kplt = seaborn.kdeplot(methpct[group], label=group, color=color)
        kplt.set(xlim=(0,1))
    plt.legend()
    plt.xlabel("% methylation")
    plt.ylabel("Density")
    plt.savefig(methpng, bbox_inches="tight", dpi=400)

if __name__ == "__main__":
    import sys
    make_graphs(sys.argv[1], sys.argv[2], sys.argv[3])
else:
    make_graphs(snakemake.input.cgfile, snakemake.output.depthpng, snakemake.output.methpng)
