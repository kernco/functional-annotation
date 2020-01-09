import os
import json

txtout = open(snakemake.output.txt, 'w')
csvout = open(snakemake.output.csv, 'w')

headers = ["Tissue", "Rep", "Depth", "EFL", "NRF", "PBC1", "PBC2", "NSC", "RSC", "JSD", "SJSD", "IJSD"]
txtout.write("{: >15s}{: >7s}{: >15s}{: >10s}{: >10s}{: >10s}{: >10s}{: >10s}{: >10s}{: >10s}{: >10s}{: >10s}\n".format(*headers))
csvout.write(','.join(headers) + '\n')

for sample in snakemake.params.libraries:
    tissue = sample.split('_')[1]
    try:
        replicate = sample.split('_')[2]
    except IndexError:
        replicate = ' '
    output = "{: >15s}{: >7s}".format(tissue, replicate)
    data = [tissue, replicate]
    with open('Metrics/{}_Alignment_Stats.json'.format(sample)) as f:
        stats = json.loads(f.read())
    output += "{: >15,d}{: >10d}{: >10.2f}{: >10.2f}{: >10.2f}".format(int(stats['Final Reads']), 0, float(stats['NRF']), float(stats['PBC1']), float(stats['PBC2']))
    data += [int(stats['Final Reads']), 0, float(stats['NRF']), float(stats['PBC1']), float(stats['PBC2'])]
    if os.path.isfile('Aligned_Reads/{}.bam'.format(sample)):
        try:
            output += "{: >10.2f}{: >10.2f}".format(float(stats['NSC']), float(stats['RSC']))
            data += [float(stats['NSC']), float(stats['RSC'])]
        except KeyError:
            output += "{: >10.2f}{: >10.2f}".format(0, 0)
            data += [0, 0]
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

