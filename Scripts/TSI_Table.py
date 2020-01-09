import math

def calc_index(nums):
    expm = max(nums)
    if expm == 0:
        return 0
    index = sum([1 - (x/expm) for x in nums]) / (len(nums) - 1)
    return index

def make_tsi_table(countfile, outfile, skip=[]):
    data = []
    with open(countfile) as f:
        headers = f.readline().split()
        keep = [x not in skip for x in headers]
        keptheaders = [x for x, y in zip(headers, keep) if y]
        for line in f:
            cols = line.split()
            #tpms = [math.log2(float(x)+1) if float(x) > 1 else 0 for x, y in zip(cols[1:], keep) if y]
            #tpms = [math.log2(float(x)+1) for x, y in zip(cols[1:], keep) if y]
            tpms = [float(x) for x, y in zip(cols[1:], keep) if y]
            tsi = calc_index(tpms)
            maxtpm = max(tpms)
            tissue = keptheaders[tpms.index(maxtpm)].split('_')[0]
            data.append([tsi, cols[0], maxtpm, tissue])

    with open(outfile, 'w') as f:
        print("\t".join(["TSI", "GeneID", "TPM", "Tissue"]), file=f)
        for line in reversed(sorted(data)):
            print('\t'.join([str(x) for x in line]), file=f)

#if __name__ == "__main__":
    #import sys
    #make_tsi_table(sys.argv[1], sys.argv[2], sys.argv[3:])
#else:
make_tsi_table(snakemake.input.table, snakemake.output.table)
