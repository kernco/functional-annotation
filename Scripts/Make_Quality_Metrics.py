import subprocess
import os

txtout = open(snakemake.output.txt, 'w')
csvout = open(snakemake.output.csv, 'w')

headers = ["Tissue", "Rep", "Depth", "NRF", "PBC1", "PBC2", "NSC", "RSC", "JSD", "SJSD", "IJSD"]
txtout.write("{: >15s}{: >7s}{: >15s}{: >10s}{: >10s}{: >10s}{: >10s}{: >10s}{: >10s}{: >10s}{: >10s}\n".format(*headers))
csvout.write(','.join(headers) + '\n')

for sample in snakemake.params.libraries:
    tissue = sample.split('_')[1]
    try:
        replicate = sample.split('_')[2]
    except IndexError:
        replicate = ' '
    output = "{: >15s}{: >7s}".format(tissue, replicate)
    data = [tissue, replicate]
    depth = int(subprocess.check_output(["samtools", "view", "-c", "Aligned_Reads/{}.bam".format(sample)]))
    result = subprocess.check_output("bamToBed -i Bwa_Output/{}.filtered.bam | awk 'BEGIN{{OFS=\"\t\"}}{{print $1,$2,$3,$6}}' | sort | uniq -c | awk 'BEGIN{{mt=1;m0=1;m1=1;m2=1}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} END{{printf \"%d\t%d\t%d\t%d\t%f\t%f\t%f\",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}}'".format(sample), shell=True)
    metrics = result.split()
    output += "{: >15,d}{: >10.2f}{: >10.2f}{: >10.2f}".format(depth, float(metrics[4]), float(metrics[5]), float(metrics[6]))
    data += [depth, float(metrics[4]), float(metrics[5]), float(metrics[6])]
    if os.path.isfile('Aligned_Reads/{}.bam'.format(sample)):
        with open('Metrics/{}_spp_stats.txt'.format(sample)) as f:
            sppstats = f.read().split()
        output += "{: >10.2f}{: >10.2f}".format(float(sppstats[8]), float(sppstats[9]))
        data += [float(sppstats[8]), float(sppstats[9])]
        with open('Metrics/{}_DeepTools_Metrics.txt'.format(sample)) as f:
            lines = f.readlines()
            chipstats = lines[1].split()
            try:
                sjsd = float(chipstats[8])
            except IndexError:
                sjsd = 0
            if len(lines) > 2:
                jsd = float(chipstats[7])
                ijsd = float(lines[2].split()[8])
                if ijsd > sjsd:
                    jsd = -jsd #Input is less uniform than ChIP
                output += "{: >10.2f}{: >10.2f}{: >10.2f}".format(jsd, sjsd, ijsd)
                data += [jsd, sjsd, ijsd]
            else:
                output += "{: >10s}{: >10.2f}{: >10s}".format(' ', sjsd, ' ')
                data += [' ', sjsd, ' ']
    else:
        output += "{: ^10s}{: ^10s}{: ^10s}".format(' ', ' ', ' ')
        data += [' ', ' ', ' ']
    txtout.write(output + "\n")
    csvout.write(','.join([str(x) for x in data]) + '\n')

