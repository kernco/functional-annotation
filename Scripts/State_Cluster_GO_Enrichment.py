from collections import defaultdict
from goatools.obo_parser import GODag
from goatools.associations import read_associations
from goatools.go_enrichment import GOEnrichmentStudy

genes = defaultdict(set)
background = set()
if snakemake.wildcards.state_type == 'Promoter':
    min_dist = 0
    max_dist = 1000
elif snakemake.wildcards.state_type == 'Enhancer':
    min_dist = 5000
    max_dist = 50000
else:
    sys.exit(-1)
with open(snakemake.input.clusters) as f:
    for line in f:
        cols = line.strip().split()
        if int(cols[-1]) <= max_dist and int(cols[-1]) >= min_dist:
            genes[cols[3]].add(cols[7])
            background.add(cols[7])

obodag = GODag("go-basic.obo")
id2go = read_associations("sym2go.txt")
goeaobj = GOEnrichmentStudy(background, id2go, obodag, propagate_counts=False, alpha=0.05, methods=['fdr_bh'])
outfile = open(snakemake.output.txt, 'w')
for cluster, geneids in sorted(genes.items()):
    outfile.write("Cluster {}\n".format(cluster))
    goea_results_all = goeaobj.run_study(geneids)
    for fdr, name, enrichment in sorted([(r.p_fdr_bh, r.name, r.enrichment) for r in goea_results_all if r.p_fdr_bh < 0.2]):
        outfile.write("\t{}\t{}\t{}\n".format(name, fdr, enrichment))
    outfile.write("\n")
    #GOEnrichmentStudy.print_summary(goea_results_sig)






