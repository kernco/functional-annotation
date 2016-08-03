import sys

peaks = set()
with open(sys.argv[1]) as f:
    for line in f:
        peaks.add(line.split()[3].split('/')[-1])

with open(sys.argv[2]) as f:
    for line in f:
        if line.split()[3].split('/')[-1] in peaks:
            print line.strip()

