import sys

tissue = sys.argv[1].split('_')[1]
for line in sys.stdin:
    if line.startswith('@SQ'):
        parts = line.split(':')
        parts[1] = tissue + "_" + parts[1]
        sys.stdout.write(':'.join(parts))
    elif line.startswith('@'):
        sys.stdout.write(line)
    else:
        parts = line.split('\t')
        parts[2] = tissue + "_" + parts[2]
        sys.stdout.write('\t'.join(parts))
