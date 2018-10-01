
outfile = open(snakemake.output.trackdb, 'w')
for fe in snakemake.input.fes:
    filename = fe.split('/')[-1]
    name = ' '.join(filename.split('_')[:-1])
    outfile.write('track ' + filename + '\n')
    outfile.write('bigDataUrl ' + filename + '\n')
    outfile.write('shortLabel ' + name + '\n')
    outfile.write('longLabel ' + name + '\n')
    outfile.write('maxHeightPixels 128:32:11\n')
    outfile.write('visibility 0\n')
    outfile.write('autoScale on\n')
    outfile.write('type bigWig\n\n')

for peak in snakemake.input.narrow:
    filename = peak.split('/')[-1]
    name = ' '.join(filename.split('_')[:-1]) + " Peaks"
    outfile.write('track ' + filename + '\n')
    outfile.write('bigDataUrl ' + filename + '\n')
    outfile.write('shortLabel ' + name + '\n')
    outfile.write('longLabel ' + name + '\n')
    outfile.write('maxHeightPixels 128:16:11\n')
    outfile.write('visibility 0\n')
    outfile.write('type bigNarrowPeak\n\n')

for peak in snakemake.input.broad:
    filename = peak.split('/')[-1]
    name = ' '.join(filename.split('_')[:-1]) + " Peaks"
    outfile.write('track ' + filename + '\n')
    outfile.write('bigDataUrl ' + filename + '\n')
    outfile.write('shortLabel ' + name + '\n')
    outfile.write('longLabel ' + name + '\n')
    outfile.write('maxHeightPixels 128:16:11\n')
    outfile.write('visibility 0\n')
    outfile.write('type bigBed\n\n')

#for seg in snakemake.input.segs:
    #filename = seg.split('/')[-1]
    #name = filename.split('_')[0] + " Chromatin States"
    #outfile.write('track ' + filename + '\n')
    #outfile.write('bigDataUrl ' + filename + '\n')
    #outfile.write('shortLabel ' + name + '\n')
    #outfile.write('longLabel ' + name + '\n')
    #outfile.write('maxHeightPixels 128:16:11\n')
    #outfile.write('visibility 0\n')
    #outfile.write('itemRgb on\n')
    #outfile.write('type bigBed\n\n')
