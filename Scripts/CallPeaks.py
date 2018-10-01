import gzip
import subprocess

def peak_type(library):
    """ Given a library or assay, returns whether the assay is 'broad' or 'narrow' """
    if "_" in library:
        assay, tissue, replicate = library.split("_")
    else:
        assay = library
    if assay in snakemake.config['broad_peaks']:
        return 'broad'
    elif assay in snakemake.config['narrow_peaks']:
        return 'narrow'
    else:
        raise Exception

def peak_file(wildcards):
    return 'Macs2/{library}_peaks.{type}Peak'.format(library=wildcards.library, type=peak_type(wildcards.library))

with gzip.open(snakemake.input.chip) as f:
    chip_depth = 0
    for line in f:
        chip_depth += 1
with gzip.open(snakemake.input.control) as f:
    input_depth = 0
    for line in f:
        input_depth += 1
if chip_depth > input_depth:
    snakemake.params.scaling = '--to-large'
else:
    snakemake.params.scaling = ''

peaktype = peak_type(snakemake.wildcards.library)
command = ['macs2', 'callpeak', '-t', snakemake.input.chip]
if hasattr(snakemake.input, 'control'):
    command += ['-c', snakemake.input.control]
command += ['-f', 'BED', '-n', 'Macs2/{}'.format(snakemake.wildcards.library), '-g', snakemake.config["genomesize"], '-B', '--SPMR', '--keep-dup', 'all', '--fix-bimodal', '--extsize', '200']
if snakemake.params.scaling:
    command.append(snakemake.params.scaling)
if peaktype == 'narrow':
    command += ['-q', '0.01']
elif peaktype == 'broad':
    command += ['-q', '0.05', '--broad']
print(command)
subprocess.run([str(arg) for arg in command])

peakfile = 'Macs2/{}_peaks.{}Peak'.format(snakemake.wildcards.library, peaktype)
subprocess.run(['cp', peakfile, snakemake.output.peaks])
