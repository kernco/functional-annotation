#Use this command on a bam file to pipe data to this script:
#samtools depth bam | cut -f3 | this script
import sys

exact1 = 0
atleast1 = 0
atleast2 = 0
for line in sys.stdin:
    if int(line) == 1:
        exact1 += 1
        atleast1 += 1
    elif int(line) == 2:
        atleast1 += 1
        atleast2 += 1
    else:
        atleast1 += 1
        atleast2 += 1

print "PBC1:", float(exact1) / atleast1
print "PBC2:", float(exact1) / atleast2
