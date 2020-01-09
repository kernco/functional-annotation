from collections import defaultdict
import statistics

#outfile = open(snakemake.output.csv, 'w')
def load_counts(expression_files):
    counts = defaultdict(lambda: defaultdict(float))
    for filename in expression_files:
        with open(filename) as f:
            headers = f.readline().split()
            rep = headers[-1].split('_')[-1]
            tissues = [x.split('_')[0] for x in headers]
            for line in f:
                cols = line.split()
                tpms = [float(x) for x in cols[1:]]
                for tpm, tissue in zip(tpms,tissues):
                    counts[tissue][cols[0]] += tpm

    for tissue, genes in counts.items():
        for gene, count in genes.items():
            counts[tissue][gene] /= len(expression_files)

    return counts

def load_annotation(gtf):
    mrna = set()
    ncrna = set()
    with open(gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.split()
            if cols[2] == 'gene':
                geneid = cols[9].strip(';').strip('"')
                if 'protein_coding' in line:
                    mrna.add(geneid)
                else:
                    ncrna.add(geneid)
    return mrna, ncrna

def print_groups(data):
    total = 0
    counter = 0
    for count, group in reversed(sorted([(y, x) for x, y in data.items()])):
        if ',' in group and counter < 10:
            print(count, group)
            counter += 1
            total += count
    print(total)

def make_report(expression_files, gtf):
    mrna, ncrna = load_annotation(gtf)
    nx = load_counts(expression_files)
    genelist = set()
    print(len(mrna), len(ncrna))
    for tissue, genes in nx.items():
        print(tissue, len([x for gene, x in genes.items() if x >= 1 and gene in mrna]), len([x for gene, x in genes.items() if x >= 1 and gene in ncrna]))
        genelist.update(genes.keys())

    mrna_detection = defaultdict(int)
    ncrna_detection = defaultdict(int)
    mrna_enriched = defaultdict(int)
    ncrna_enriched = defaultdict(int)
    mrna_enhanced = defaultdict(int)
    ncrna_enhanced = defaultdict(int)
    testset = []
    for gene in genelist:
        expression = list(reversed(sorted([(nx[tissue][gene], tissue) for tissue in nx])))
        num_detected = len([x for x, t in expression if x >= 1])
        if gene in mrna:
            mrna_detection[num_detected] += 1
        elif gene in ncrna:
            ncrna_detection[num_detected] += 1
        if num_detected > 0:
            for i in range(len(expression)-3):
                if expression[i][0] > expression[i+1][0]*4:
                    if gene in mrna:
                        mrna_enriched[','.join(sorted([x[1] for x in expression[:i+1]]))] += 1
                        if i == 0 and expression[0][1] == 'Adipose':
                            testset.append(gene)
                    elif gene in ncrna:
                        ncrna_enriched[','.join(sorted([x[1] for x in expression[:i+1]]))] += 1
                    break
            else:
                if expression[0][0] > statistics.mean(exp[0] for exp in expression[1:])*4:
                    if gene in mrna:
                        mrna_enhanced[expression[0][1]] += 1
                    elif gene in ncrna:
                        ncrna_enhanced[expression[0][1]] += 1

    print(mrna_detection)
    print(ncrna_detection)
    for tissue in nx.keys():
        print(tissue, mrna_enriched[tissue], ncrna_enriched[tissue], mrna_enhanced[tissue], ncrna_enhanced[tissue])
    print_groups(mrna_enriched)
    print_groups(ncrna_enriched)
    for gene in testset:
        print(gene)


    #print("Detected genes in " + filename, file=outfile)
    #print(','.join(["Tissue"] + [str(cutoff) for cutoff in cutoffs]), file=outfile)
    #for tissue in tissues:
        #print(','.join([tissue] + [str(counts[tissue][cutoff]) for cutoff in cutoffs]), file=outfile)
    #print("", file=outfile)
#
#print("Detected in all replicates", file=outfile)
#print(','.join(["Tissue"] + [str(cutoff) for cutoff in cutoffs]), file=outfile)
#for tissue in tissues:
    #nums = []
    #for cutoff in cutoffs:
        #genes = None
        #for filename in snakemake.input.tpms:
            #if not genes:
                #genes = expressed[filename][tissue][cutoff]
            #else:
                #genes &= expressed[filename][tissue][cutoff]
        #nums.append(str(len(genes)))
    #print(','.join([tissue] + nums), file=outfile)
#
#print("", file=outfile)

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

#if __name__ == "__main__":
   #import sys
   #make_report(sys.argv[2:], sys.argv[1])
#else:
make_report(snakemake.input.tpms, snakemake.input.annotation)
