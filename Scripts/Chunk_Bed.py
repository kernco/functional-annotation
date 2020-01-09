
size = 200
outfile = open(snakemake.output.outfile, 'w')
with open(snakemake.input.infile) as f:
    for line in f:
        cols = line.strip().split()
        start = int(cols[1])
        nums = list(range(int(cols[1]), int(cols[2]), size)) + [int(cols[2])]
        for i,j in zip(nums, nums[1:]):
            cols[1] = str(i)
            cols[2] = str(j)
            outfile.write("\t".join(cols) + "\n")
