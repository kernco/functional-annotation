import matplotlib.pyplot as plt

data = []
for filename in reversed(snakemake.input):
    with open(filename) as f:
        data.append([float(line.split()[-1].strip(';').strip('"')) for line in f])

fig = plt.figure()
plt.boxplot(data, sym='', vert=False, labels=list(reversed(snakemake.params.labels)))
plt.xlabel('Expression (TPM)')
fig.savefig(snakemake.output.png, bbox_inches='tight')
plt.close(fig)

