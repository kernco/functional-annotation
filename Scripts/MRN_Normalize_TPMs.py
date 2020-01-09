from collections import defaultdict
import numpy
import statistics
import math

tpm = defaultdict(dict)
samples = set()
groups = set()
replicates = set()
symbols = {}

def load_data(filenames):
    for filename in filenames:
        with open(filename) as f:
            headers = f.readline().split()
            for sample in headers[2:]:
                samples.add(sample)
                group, replicate = sample.split('_')
                groups.add(group)
                replicates.add(replicate)
            for line in f:
                cols = line.split()
                symbols[cols[0]] = cols[1]
                for val, header in zip(cols[2:], headers[2:]):
                    tpm[cols[0]][header] = float(val)

def mean_exp_in_group(gene, group):
    vals = [val for sample, val in tpm[gene].items() if sample.split('_')[0] == group]
    return sum(vals) / len(vals)

#Load TPMs from tables
load_data(snakemake.input.tpmtables)

#Remove genes with TPM<=cutoff in all libraries
cutoff = 0
genes = [gene for gene, tpms in tpm.items() if max(tpms.values()) <= cutoff]
for gene in genes:
    del tpm[gene]

#Pre-normalize by library size
N = {}
for sample in samples:
    N[sample] = sum([tpm[gene][sample] for gene in tpm])

for gene, tpms in tpm.items():
    for sample in tpms:
        tpm[gene][sample] /= N[sample]

#Find representative group
quartiles = defaultdict(list)
for sample in samples:
    quartiles[sample.split('_')[0]].append(numpy.percentile([tpm[gene][sample] for gene in tpm], 75))
mean = sum([sum(qts) / len(qts) for sample, qts in quartiles.items()]) / len(quartiles)
rep_group = sorted([(abs(mean - (sum(qts) / len(qts))), sample) for sample, qts in quartiles.items()])[0][1]
rep_group = rep_group.split('_')[0]

#Representative expression
Y_gene = {}
for gene, vals in tpm.items():
    Y_gene[gene] = mean_exp_in_group(gene, rep_group)

#Relative sizes of transcriptomes and representative group
tau = {}
for group in groups:
    tau[group] = statistics.median([mean_exp_in_group(gene, group) / Y_gene[gene] for gene in tpm if Y_gene[gene] > 0])

#Effective library size
e = {}
for sample in samples:
    e[sample] = tau[sample.split('_')[0]] * N[sample]

#Normalization factors
f = {}
ft = math.exp(sum([math.log(x) for x in e.values()]) / len(samples))
for sample in samples:
    f[sample] = e[sample] / ft

#Normalize TPMs
for gene, tpms in tpm.items():
    for sample, val in tpms.items():
        tpm[gene][sample] = (val / f[sample]) * 1000000

for filename in snakemake.output.tpmmrns:
    replicate = filename.split('.')[0].split('_')[-1]
    subsamples = [sample for sample in sorted(samples) if sample.split('_')[-1] == replicate]
    with open(filename, "w") as outfile:
        print("\t".join(["Gene", "RefID"] + sorted(subsamples)), file=outfile)
        for gene, tpms in tpm.items():
            nums = [str(tpms[sample]) for sample in subsamples]
            print("\t".join([gene, symbols[gene]] + nums), file=outfile)

