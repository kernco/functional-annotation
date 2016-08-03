import sys

with open(sys.argv[1]) as f:
    name = None
    seq = ''
    flag = False
    for line in f:
        if line.startswith(">"):
            if line.strip().split('.')[-1] == 'exon1':
                if flag and len(seq) > 200:
                    sys.stdout.write(name + "\n" + seq + "\n")
                name = '.'.join(line.split('.')[:-1])
                flag = False
                seq = ''
                #print '\n' + '.'.join(line.split('.')[:-1])
            else:
                flag = True
        else:
            seq += line.strip()
            #sys.stdout.write(line.strip())

if flag and len(seq) > 200:
    sys.stdout.write(name + "\n" + seq)
