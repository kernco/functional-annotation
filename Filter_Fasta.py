import sys

remove = set()
with open(sys.argv[1]) as f:
    for line in f:
        remove.add(line.split()[0].split('.')[0])

with open(sys.argv[2]) as f:
    flag = False
    for line in f:
        if line.startswith('>'):
            if line.strip()[1:].split('.')[0] not in remove:
                flag = True
            else:
                flag = False
        if flag:
            sys.stdout.write(line)


