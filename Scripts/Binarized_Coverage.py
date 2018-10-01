import sys
from collections import defaultdict

counts = defaultdict(int)
total = 0
with open(sys.argv[1]) as f:
    f.readline()
    headers = f.readline().split()
    for line in f:
        nums = line.split()
        for i, num in enumerate(nums):
            if num == '1':
                counts[i] += 1
        total += 1

for k, v in counts.items():
    print(headers[k], v / total)

