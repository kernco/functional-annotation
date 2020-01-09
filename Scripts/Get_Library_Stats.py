import json
import subprocess

with open(snakemake.input.json) as f:
    stats = json.loads(f.read())


stats['Peak Calls'] = int(subprocess.check_output(["wc", "-l", snakemake.input.peaks]).split()[0])
with open(snakemake.input.frip) as f:
    f.readline()
    cols = f.readline().split()
    stats['FRiP'] = float(cols[3]) / float(cols[4])
    stats['Signal Reads'] = int(stats['FRiP'] * stats['Final Reads'])
    stats['Signal Reads per KB of Genome'] = (stats['Signal Reads'] * 1000) / snakemake.config['genomesize']

with open(snakemake.output.json, 'w') as f:
    f.write(json.dumps(stats, indent=4, separators=(',', ': ')) + '\n')
