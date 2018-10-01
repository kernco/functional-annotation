
######## Snakemake header ########
import sys; sys.path.append("/share/apps/bio3user/miniconda3/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X\'\x00\x00\x00Track_Hub/H3K4me3_Spleen_M08_Sorted.bdgq\x06X\x1f\x00\x00\x00ChromHMM/Chromosome_Lengths.txtq\x07e}q\x08(X\x06\x00\x00\x00_namesq\t}q\n(X\x03\x00\x00\x00bedq\x0bK\x00N\x86q\x0cX\n\x00\x00\x00chromsizesq\rK\x01N\x86q\x0euh\x0bh\x06h\rh\x07ubX\x06\x00\x00\x00outputq\x0fcsnakemake.io\nOutputFiles\nq\x10)\x81q\x11X/\x00\x00\x00Track_Hub/H3K4me3_Spleen_M08_Sorted.bdg_trimmedq\x12a}q\x13(h\t}q\x14X\x07\x00\x00\x00trimmedq\x15K\x00N\x86q\x16sh\x15h\x12ubX\x06\x00\x00\x00paramsq\x17csnakemake.io\nParams\nq\x18)\x81q\x19}q\x1ah\t}q\x1bsbX\t\x00\x00\x00wildcardsq\x1ccsnakemake.io\nWildcards\nq\x1d)\x81q\x1eX\'\x00\x00\x00Track_Hub/H3K4me3_Spleen_M08_Sorted.bdgq\x1fa}q (h\t}q!X\x06\x00\x00\x00prefixq"K\x00N\x86q#sX\x06\x00\x00\x00prefixq$h\x1fubX\x07\x00\x00\x00threadsq%K\x01X\t\x00\x00\x00resourcesq&csnakemake.io\nResources\nq\')\x81q((K\x01K\x01e}q)(h\t}q*(X\x06\x00\x00\x00_coresq+K\x00N\x86q,X\x06\x00\x00\x00_nodesq-K\x01N\x86q.uh+K\x01h-K\x01ubX\x03\x00\x00\x00logq/csnakemake.io\nLog\nq0)\x81q1}q2h\t}q3sbX\x06\x00\x00\x00configq4}q5(X\x07\x00\x00\x00tissuesq6]q7(X\x05\x00\x00\x00Liverq8X\x04\x00\x00\x00Lungq9X\x06\x00\x00\x00Spleenq:X\n\x00\x00\x00Cerebellumq;X\x0c\x00\x00\x00Hypothalamusq<eX\x04\x00\x00\x00repsq=]q>(X\x03\x00\x00\x00M08q?X\x03\x00\x00\x00M22q@eX\x0c\x00\x00\x00narrow_peaksqA]qB(X\x07\x00\x00\x00H3K4me3qCX\x07\x00\x00\x00H3K27acqDeX\x0b\x00\x00\x00broad_peaksqE]qF(X\x08\x00\x00\x00H3K27me3qGX\x07\x00\x00\x00H3K4me1qHX\x04\x00\x00\x00ATACqIeX\x08\x00\x00\x00no_inputqJ]qK(X\x04\x00\x00\x00ATACqLX\x06\x00\x00\x00RNASeqqMeX\x06\x00\x00\x00inputsqN]qO(X\x05\x00\x00\x00InputqPX\x06\x00\x00\x00Input2qQeX\r\x00\x00\x00H3K4me3_inputqRX\x05\x00\x00\x00InputqSX\x0e\x00\x00\x00H3K27me3_inputqTX\x05\x00\x00\x00InputqUX\r\x00\x00\x00H3K4me1_inputqVX\x05\x00\x00\x00InputqWX\r\x00\x00\x00H3K27ac_inputqXX\x05\x00\x00\x00InputqYX\x0e\x00\x00\x00override_inputqZccollections\nOrderedDict\nq[)Rq\\(X\x11\x00\x00\x00H3K4me3_Liver_M08q]X\x10\x00\x00\x00Input2_Liver_M08q^X\x12\x00\x00\x00H3K27me3_Liver_M08q_X\x10\x00\x00\x00Input2_Liver_M08q`X\x11\x00\x00\x00H3K4me3_Liver_M22qaX\x10\x00\x00\x00Input2_Liver_M22qbX\x12\x00\x00\x00H3K27me3_Liver_M22qcX\x10\x00\x00\x00Input2_Liver_M22qdX\x11\x00\x00\x00H3K27me3_Lung_M22qeX\x0f\x00\x00\x00Input2_Lung_M22qfX\x13\x00\x00\x00H3K27me3_Spleen_M08qgX\x11\x00\x00\x00Input2_Spleen_M08qhX\x12\x00\x00\x00H3K4me1_Spleen_M08qiX\x11\x00\x00\x00Input2_Spleen_M08qjX\x12\x00\x00\x00H3K4me3_Spleen_M08qkX\x11\x00\x00\x00Input2_Spleen_M08qlX\x13\x00\x00\x00H3K27me3_Spleen_M22qmX\x11\x00\x00\x00Input2_Spleen_M22qnX\x12\x00\x00\x00H3K4me1_Spleen_M22qoX\x11\x00\x00\x00Input2_Spleen_M22qpuX\x06\x00\x00\x00genomeqqX)\x00\x00\x00/group/zhougrp/Genomes/bosTau8/bosTau8.faqrX\n\x00\x00\x00genomesizeqs\x8a\x05\x13\xf4\xef\x9d\x00X\n\x00\x00\x00annotationqtX4\x00\x00\x00/group/zhougrp/Genomes/bosTau8/bosTau8_Ensembl91.gtfquX\x04\x00\x00\x00mapqqvK\x1eX\n\x00\x00\x00chromsizesqwX\x1f\x00\x00\x00ChromHMM/Chromosome_Lengths.txtqxX\x0f\x00\x00\x00ChromHMM_genomeqyX\x07\x00\x00\x00bosTau8qzuX\x04\x00\x00\x00ruleq{X\x08\x00\x00\x00Trim_Bedq|ub.'); from snakemake.logging import logger; logger.printshellcmds = True
######## Original script #########

lens = {}
with open(snakemake.input.chromsizes) as f:
    for line in f:
        cols = line.split()
        lens[cols[0]] = int(cols[1])

outfile = open(snakemake.output.trimmed, "w")
with open(snakemake.input.bed) as f:
    for line in f:
        cols = line.split()
        if int(cols[2]) > lens[cols[0]]:
            cols[2] = str(lens[cols[0]])
        outfile.write('\t'.join(cols) + "\n")
