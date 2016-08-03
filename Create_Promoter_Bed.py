import sys

promoters = {}
with open(sys.argv[1]) as f:
    for line in f:
        cols = line.strip().split('\t')
        if cols[5] == '+':
            cols[1] = str(max(0, int(cols[1]) - 2000))
            cols[2] = cols[1]
        elif cols[5] == '-':
            cols[1] = cols[2]
            cols[2] = str(int(cols[2]) + 2000)
        print '\t'.join(cols)




