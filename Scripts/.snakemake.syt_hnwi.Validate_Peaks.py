
######## Snakemake header ########
import sys; sys.path.append("/share/apps/bio3user/miniconda3/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05X=\x00\x00\x00Macs2/H3K4me3_Cerebellum_B_Peak_Regions_With_Replicate_FE.txtq\x06a}q\x07(X\x06\x00\x00\x00_namesq\x08}q\tX\x05\x00\x00\x00peaksq\nK\x00N\x86q\x0bsh\nh\x06ubX\x06\x00\x00\x00outputq\x0ccsnakemake.io\nOutputFiles\nq\r)\x81q\x0eX;\x00\x00\x00Macs2/H3K4me3_Cerebellum_B_Peaks_Validated_by_Replicate.bedq\x0fa}q\x10(h\x08}q\x11X\x07\x00\x00\x00outfileq\x12K\x00N\x86q\x13sh\x12h\x0fubX\x06\x00\x00\x00paramsq\x14csnakemake.io\nParams\nq\x15)\x81q\x16X\x06\x00\x00\x00narrowq\x17a}q\x18(h\x08}q\x19X\x08\x00\x00\x00peaktypeq\x1aK\x00N\x86q\x1bsh\x1ah\x17ubX\t\x00\x00\x00wildcardsq\x1ccsnakemake.io\nWildcards\nq\x1d)\x81q\x1eX\x14\x00\x00\x00H3K4me3_Cerebellum_Bq\x1fa}q (h\x08}q!X\x07\x00\x00\x00libraryq"K\x00N\x86q#sX\x07\x00\x00\x00libraryq$h\x1fubX\x07\x00\x00\x00threadsq%K\x01X\t\x00\x00\x00resourcesq&csnakemake.io\nResources\nq\')\x81q((K\x01K\x01e}q)(h\x08}q*(X\x06\x00\x00\x00_coresq+K\x00N\x86q,X\x06\x00\x00\x00_nodesq-K\x01N\x86q.uh+K\x01h-K\x01ubX\x03\x00\x00\x00logq/csnakemake.io\nLog\nq0)\x81q1}q2h\x08}q3sbX\x06\x00\x00\x00configq4}q5(X\x07\x00\x00\x00tissuesq6]q7(X\x05\x00\x00\x00Liverq8X\x04\x00\x00\x00Lungq9X\x06\x00\x00\x00Spleenq:X\n\x00\x00\x00Cerebellumq;eX\x04\x00\x00\x00repsq<]q=(X\x01\x00\x00\x00Aq>X\x01\x00\x00\x00Bq?eX\x0c\x00\x00\x00narrow_peaksq@]qA(X\x07\x00\x00\x00H3K4me3qBX\x07\x00\x00\x00H3K27acqCX\x04\x00\x00\x00CTCFqDeX\x0b\x00\x00\x00broad_peaksqE]qF(X\x08\x00\x00\x00H3K27me3qGX\x07\x00\x00\x00H3K4me1qHX\x08\x00\x00\x00DNaseSeqqIeX\x08\x00\x00\x00no_inputqJ]qK(X\x08\x00\x00\x00DNaseSeqqLX\x06\x00\x00\x00RNASeqqMeX\x06\x00\x00\x00inputsqN]qO(X\x05\x00\x00\x00InputqPX\x06\x00\x00\x00Input2qQX\t\x00\x00\x00InputCTCFqRX\x08\x00\x00\x00NewInputqSeX\x0e\x00\x00\x00override_inputqTccollections\nOrderedDict\nqU)RqV(X\x10\x00\x00\x00H3K27ac_Spleen_AqWX\x0e\x00\x00\x00Input_Spleen_AqXX\x10\x00\x00\x00H3K4me1_Spleen_AqYX\x0e\x00\x00\x00Input_Spleen_AqZX\x10\x00\x00\x00H3K4me1_Spleen_Bq[X\x0e\x00\x00\x00Input_Spleen_Bq\\X\x14\x00\x00\x00H3K27ac_Cerebellum_Aq]X\x12\x00\x00\x00Input_Cerebellum_Aq^X\x14\x00\x00\x00H3K27ac_Cerebellum_Bq_X\x12\x00\x00\x00Input_Cerebellum_Bq`X\x14\x00\x00\x00H3K4me1_Cerebellum_AqaX\x12\x00\x00\x00Input_Cerebellum_AqbX\x14\x00\x00\x00H3K4me1_Cerebellum_BqcX\x12\x00\x00\x00Input_Cerebellum_BqdX\x10\x00\x00\x00H3K27ac_Cortex_BqeX\x0e\x00\x00\x00Input_Cortex_BqfX\x10\x00\x00\x00H3K4me1_Cortex_BqgX\x0e\x00\x00\x00Input_Cortex_BqhuX\r\x00\x00\x00H3K4me3_inputqiX\x05\x00\x00\x00InputqjX\x0e\x00\x00\x00H3K27me3_inputqkX\x05\x00\x00\x00InputqlX\r\x00\x00\x00H3K4me1_inputqmX\x06\x00\x00\x00Input2qnX\r\x00\x00\x00H3K27ac_inputqoX\x06\x00\x00\x00Input2qpX\n\x00\x00\x00CTCF_inputqqX\t\x00\x00\x00InputCTCFqrX\x06\x00\x00\x00genomeqsX)\x00\x00\x00/group/zhougrp/Genomes/galGal5/galGal5.faqtX\n\x00\x00\x00genomesizequJ\r\xdc\xa0HX\n\x00\x00\x00annotationqvX4\x00\x00\x00/group/zhougrp/Genomes/galGal5/galGal5_Ensembl91.gtfqwX\x04\x00\x00\x00mapqqxK\x1eX\x0f\x00\x00\x00ChromHMM_genomeqyX\x07\x00\x00\x00galGal5qzX\n\x00\x00\x00chromsizesq{X\x1f\x00\x00\x00ChromHMM/Chromosome_Lengths.txtq|uX\x04\x00\x00\x00ruleq}X\x0e\x00\x00\x00validate_peaksq~ub.'); from snakemake.logging import logger; logger.printshellcmds = True
######## Original script #########
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

