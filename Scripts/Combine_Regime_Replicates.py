from collections import defaultdict

locs = defaultdict(list)

def loadfile(filename, tag):
    with open(filename) as f:
        for line in f:
            cols = line.split()
            loc = (cols[0], int(cols[1]), int(cols[2]))#'\t'.join(cols[0:3])
            locs[loc].append(tag)

loadfile(snakemake.input.regime1[0], ("Rep1", "Reg1"))
loadfile(snakemake.input.regime2[0], ("Rep1", "Reg2"))
loadfile(snakemake.input.regime1[1], ("Rep2", "Reg1"))
loadfile(snakemake.input.regime2[1], ("Rep2", "Reg2"))

outreg1 = []# open(snakemake.output.combinereg1, "w")
outreg2 = []#open(snakemake.output.combinereg2, "w")
for loc, tags in locs.items():
    if len(tags) == 1:
        if tags[0][1] == "Reg2": #Strong vs. none = weak
            #print(loc, file=outreg1)
            outreg1.append(loc)
    else:
        if tags[0][1] == "Reg2" and tags[1][1] == "Reg2": #Strong vs. strong = strong
            #print(loc, file=outreg2)
            outreg2.append(loc)
        elif tags[0][1] == "Reg1" and tags[1][1] == "Reg1": #Weak vs. weak = weak
            #print(loc, file=outreg1)
            outreg1.append(loc)
        elif tags[0][1] == "Reg1" and tags[1][1] == "Reg2": #Strong vs. weak = weak
            #print(loc, file=outreg1)
            outreg1.append(loc)
        elif tags[0][1] == "Reg2" and tags[1][1] == "Reg1": #Strong vs. weak = weak
            #print(loc, file=outreg1)
            outreg1.append(loc)

with open(snakemake.output.combinereg1, "w") as f:
    for loc in sorted(outreg1):
        print('\t'.join([str(x) for x in loc]), file=f)

with open(snakemake.output.combinereg2, "w") as f:
    for loc in sorted(outreg2):
        print('\t'.join([str(x) for x in loc]), file=f)
