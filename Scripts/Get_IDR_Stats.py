import json
import subprocess

def lines(filename):
    return subprocess.check_output(["wc", "-l", filename]).split()[0]

stats = {}
stats["Pooled Peaks"] = int(lines(snakemake.input.pooled_peaks))
stats["True pooled"] = int(lines(snakemake.input.true_pooled))
stats["Pseudo pooled"] = int(lines(snakemake.input.pseudo_pooled))
stats["Stable peaks"] = int(lines(snakemake.input.stable_peaks))
stats["Rep1 self"] = int(lines(snakemake.input.rep1_self))
stats["Rep2 self"] = int(lines(snakemake.input.rep2_self))
Nt = stats["True pooled"]
Np = stats["Pseudo pooled"]
N1 = stats["Rep1 self"]
N2 = stats["Rep2 self"]
stats["Rescue ratio"] = max(Np, Nt) / min(Np, Nt)
stats["Self-consistency ratio"] = max(N1, N2) / min(N1, N2)

with open(snakemake.output.json, 'w') as f:
    f.write(json.dumps(stats, indent=4, separators=(',', ': ')) + '\n')
