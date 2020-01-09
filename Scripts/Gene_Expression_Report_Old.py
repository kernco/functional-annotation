from collections import defaultdict

cutoffs = [0, 0.5, 1, 2]

outfile = open(snakemake.output.csv, 'w')
expressed = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
for filename in snakemake.input.tpms:
    with open(filename) as f:
        headers = f.readline().split()
        rep = headers[-1].split('_')[-1]
        #tissues = [x.split('_')[0] for x in headers[2:]]
        tissues = [x.split('_')[0] for x in headers]
        counts = defaultdict(lambda: defaultdict(int))
        for line in f:
            cols = line.split()
            #tpms = [float(x) for x in cols[2:]]
            tpms = [float(x) for x in cols[1:]]
            for tpm, tissue in zip(tpms,tissues):
                for cutoff in cutoffs:
                    if tpm > cutoff:
                        counts[tissue][cutoff] += 1
                        expressed[filename][tissue][cutoff].add(cols[0])
    print("Detected genes in " + filename, file=outfile)
    print(','.join(["Tissue"] + [str(cutoff) for cutoff in cutoffs]), file=outfile)
    for tissue in tissues:
        print(','.join([tissue] + [str(counts[tissue][cutoff]) for cutoff in cutoffs]), file=outfile)
    print("", file=outfile)

print("Detected in all replicates", file=outfile)
print(','.join(["Tissue"] + [str(cutoff) for cutoff in cutoffs]), file=outfile)
for tissue in tissues:
    nums = []
    for cutoff in cutoffs:
        genes = None
        for filename in snakemake.input.tpms:
            if not genes:
                genes = expressed[filename][tissue][cutoff]
            else:
                genes &= expressed[filename][tissue][cutoff]
        nums.append(str(len(genes)))
    print(','.join([tissue] + nums), file=outfile)

print("", file=outfile)

#expressed = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
#for filename in snakemake.input.tsis:
    #with open(filename) as f:
        #headers = f.readline()
        #counts = defaultdict(lambda: defaultdict(int))
        #for line in f:
            #cols = line.split()
            #if float(cols[0]) >= 0.9:
                #tissue = cols[-1]
                #for cutoff in cutoffs:
                    #if float(cols[-2]) > cutoff:
                        #counts[tissue][cutoff] += 1
                        #expressed[filename][tissue][cutoff].add(cols[1])
    #print("Tissue-specific genes in " + filename, file=outfile)
    #print(','.join(["Tissue"] + [str(cutoff) for cutoff in cutoffs]), file=outfile)
    #for tissue in tissues:
        #print(','.join([tissue] + [str(counts[tissue][cutoff]) for cutoff in cutoffs]), file=outfile)
    #print("", file=outfile)

#print("Tissue-specific in all replicates", file=outfile)
#print(','.join(["Tissue"] + [str(cutoff) for cutoff in cutoffs]), file=outfile)
#for tissue in tissues:
    #nums = []
    #for cutoff in cutoffs:
        #genes = None
        #for filename in snakemake.input.tsis:
            #if not genes:
                #genes = expressed[filename][tissue][cutoff]
            #else:
                #genes &= expressed[filename][tissue][cutoff]
        #nums.append(str(len(genes)))
    #print(','.join([tissue] + nums), file=outfile)

