import subprocess

assay = snakemake.wildcards.assay
tissues = sorted(snakemake.params.tissues)
reps = sorted(snakemake.params.reps)

gsize = float(snakemake.config['genomesize'])
txtout = open(snakemake.output.txt, "w")
print("{: >12s}{: >12s}{: >9s}".format("Tissue", "Peaks", "Cov.") + ''.join(["{: >15s}".format(rep + " FRiP") for rep in reps]), file=txtout)
for tissue in tissues:
    outline = "{: >12s}".format(tissue)
    peakfile = 'Peak_Calls/{}_{}_IDR.{}Peak'.format(assay, tissue, snakemake.params.peaktype)
    fripfile = 'Metrics/{}_{}_IDR_FRiP.txt'.format(assay, tissue)
    try:
        output = subprocess.check_output(['wc', '-l', peakfile])
        total = sum((abs(int(line.split()[2]) - int(line.split()[1])) for line in open(peakfile)))
        outline += "{: >12,d}{: >8.2f}%".format(int(output.split()[0]), (total / gsize) * 100)
    except (subprocess.CalledProcessError, FileNotFoundError):
        outline += "{: >12,d}{: >8.2f}%".format(0, 0)
    for rep in reps:
        try:
            with open(fripfile) as f:
                for line in f:
                    cols = line.split()
                    if rep in cols[0]:
                        frip = float(cols[2])
            outline += "{: >14.2f}%".format(frip)
        except FileNotFoundError:
            outline += "{: >14.2f}%".format(0)
    print(outline, file=txtout)
