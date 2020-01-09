import os
from snakemake.shell import shell

if os.path.getsize(snakemake.input.peaks) > 0:
    shell('plotEnrichment --bamfiles {snakemake.input.bam} --BED {snakemake.input.peaks} -o {snakemake.output.figure} --outRawCounts {snakemake.output.metrics} -e 200 -p {snakemake.threads}')
else:
    shell('touch {snakemake.output.metrics} {snakemake.output.figure}')
