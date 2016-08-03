import sys
import re

ids = set()
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('>'):
            ids.add(line.strip().split('.')[1])

with open(sys.argv[2]) as f:
    for line in f:
        if re.search(r'transcript_id "(TCONS_[0-9]+)"', line).group(1) in ids:
            print line.strip()
