import subprocess
import json

stats = {}

stats["Assay Type"] = snakemake.params.assay_type

if snakemake.params.assay_type == 'Single':
    divfactor = 1
elif snakemake.params.assay_type in ['Paired', 'ATAC', 'RNA-seq']:
    divfactor = 2

#Get number of raw reads from TrimGalore report
with open(snakemake.input.trim_report) as f:
    for line in f:
        if line.startswith("Total reads processed"):
            stats["Raw Reads"] = int(line.split()[3].replace(',',''))

#Count number of reads in trimmed fastq
zcat = subprocess.Popen(("zcat", snakemake.input.trimmed_fq), stdout=subprocess.PIPE)
stats["Trimmed Reads"] = int(int(subprocess.check_output(('wc', '-l'), stdin=zcat.stdout)) / 4)


#Count reads in aligned, filtered, and deduplicated files
if snakemake.params.assay_type == 'RNA-seq':
    uniquely = subprocess.check_output(["grep", "Uniquely mapped reads number", snakemake.input.star_log])
    multimap = subprocess.check_output(["grep", "Number of reads mapped to multiple loci", snakemake.input.star_log])
    stats["Aligned Reads"] = int(uniquely.split(b'|')[1]) + int(multimap.split(b'|')[1])
    stats["Filtered Reads"] = int(int(subprocess.check_output(["samtools", "view", "-c", snakemake.input.final_bam]))) // 2
    stats["Final Reads"] = stats["Filtered Reads"]
else:
    stats["Aligned Reads"] = int(int(subprocess.check_output(["samtools", "view", "-c", "-F", "4", snakemake.input.aligned_bam])) / divfactor)
    stats["Filtered Reads"] = int(int(subprocess.check_output(["samtools", "view", "-c", "-F", "4", snakemake.input.filtered_bam])) / divfactor)
    stats["Deduplicated Reads"] = int(int(subprocess.check_output(["samtools", "view", "-c", snakemake.input.deduped_bam])) / divfactor)
    stats["Final Reads"] = stats["Deduplicated Reads"]
    #Get quality metrics
    result = subprocess.check_output("bamToBed -i {} | awk 'BEGIN{{OFS=\"\t\"}}{{print $1,$2,$3,$6}}' | sort | uniq -c | awk 'BEGIN{{mt=1;m0=1;m1=1;m2=1}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} END{{printf \"%d\t%d\t%d\t%d\t%f\t%f\t%f\",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}}'".format(snakemake.input.filtered_bam), shell=True)
    metrics = result.split()
    stats["NRF"] = float(metrics[4])
    stats["PBC1"] = float(metrics[5])
    stats["PBC2"] = float(metrics[6])
    try:
        with open(snakemake.input.spp_stats) as f:
            sppstats = f.read().split()
            stats['NSC'] = sppstats[8]
            stats['RSC'] = sppstats[9]
            stats['Est. Fragment Length'] = sppstats[2].split(',')[0]
    except (IndexError, AttributeError):
        stats['NSC'] = "NA"
        stats['RSC'] = "NA"
        stats['Est. Fragment Length'] = "NA"

#Write stats to file
with open(snakemake.output.json, 'w') as f:
    f.write(json.dumps(stats, indent=4, separators=(',', ': ')) + '\n')
