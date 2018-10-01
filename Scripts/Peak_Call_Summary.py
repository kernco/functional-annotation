import subprocess

assay = snakemake.wildcards.assay
tissues = snakemake.params.tissues
reps = snakemake.params.reps

gsize = float(snakemake.config['genomesize'])
txtout = open(snakemake.output.txt, "w")
csvout = open(snakemake.output.csv, "w")
print("{: >12s}".format("Tissue") + ''.join(["{: >12s} {: >7s} {:>7s}".format(rep, "Cov.", "FRiP") for rep in sorted(reps)]) + "{: >12s}".format("Combined") + "{: >12s}".format("Coverage") + ''.join(["{: >11s}".format(rep + " FRiP") for rep in sorted(reps)]), file=txtout)
csvout.write(','.join(["Tissue", "Replicate", "Peaks", "Coverage", "FRiP"]) + "\n")
for tissue in sorted(tissues):
    outline = "{: >12s}".format(tissue)
    for rep in sorted(reps):
        repid = '{}_{}_{}'.format(assay, tissue, rep)
        peakfile = 'Peak_Calls/{}_Peaks.bed'.format(repid)
        try:
            output = subprocess.check_output(['wc', '-l', peakfile])
            total = sum([abs(int(line.split()[2]) - int(line.split()[1])) for line in open(peakfile)])
            with open('Metrics/{}_{}_{}_FRiP.txt'.format(assay, tissue, rep)) as f:
                for line in f:
                    cols = line.split()
                    if rep in cols[0]:
                        frip = float(cols[2])
            outline += "{: >12,d} {: >6.2f}% {: >6.2f}%".format(int(output.split()[0]), (total / gsize) * 100, frip)
            csvout.write(','.join([tissue, str(rep), str(int(output.split()[0])), str(total / gsize), str(frip / 100)]) + '\n')
        except subprocess.CalledProcessError:
            outline += "{: >12s} {: >7s} {: >7s}".format("-", "-", "-")
    peakfile = 'Peak_Calls/{}_{}_Combined_Peaks.bed'.format(assay, tissue)
    try:
        output = subprocess.check_output(['wc', '-l', peakfile])
        outline += "{: >12,d}".format(int(output.split()[0]))
    except subprocess.CalledProcessError:
        outline += "{: >12s}".format("N/A")
    total = sum([abs(int(line.split()[2]) - int(line.split()[1])) for line in open(peakfile)])
    outline += "{: >11.2f}%".format((total / float(gsize))*100)
    print(','.join([tissue + " Combined", str(int(output.split()[0])), str(total / float(gsize))]))
    fripfile = 'Metrics/{}_{}_FRiP.txt'.format(assay, tissue)
    for rep in sorted(reps):
        with open(fripfile) as f:
            for line in f:
                cols = line.split()
                if rep in cols[0]:
                    frip = float(cols[2])
                    break
            else:
                frip = 0
        outline += "{: >10.2f}%".format(frip)
    print(outline, file=txtout)

