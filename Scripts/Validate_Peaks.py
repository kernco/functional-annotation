from collections import defaultdict

peaks = defaultdict(list)
outfile = open(snakemake.output.outfile, "w")
with open(snakemake.input.peaks) as f:
    for line in f:
        parts = line.strip().split()
        if snakemake.params.peaktype == 'broad':
            peaks[(parts[0], parts[1], parts[2], parts[3], parts[4])].append((float(parts[12]), int(parts[13])))
        else:
            peaks[(parts[0], parts[1], parts[2], parts[3], parts[4])].append((float(parts[13]), int(parts[14])))

for k, v in peaks.items():
    try:
        score = sum([x[0] * x[1] for x in v]) / sum([x[1] for x in v])
    except ZeroDivisionError:
        score = 0
    if snakemake.params.peaktype == 'broad' and score >= 1.5 or score >= 2:
        print("{}\t{}\t{}\t{}\t{}\t.".format(*k), file=outfile)

