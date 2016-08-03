import sys

seqs = {}
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith(">"):
            name = line.split()[0][1:]
            seqs[name] = 0
        else:
            seqs[name] += len(line.strip())

for tissue in sys.argv[2:]:
    for seq in seqs:
        print "{}_{}\t{}".format(tissue, seq, seqs[seq])
