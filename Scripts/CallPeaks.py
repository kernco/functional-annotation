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

def get_depth(filename):
    with gzip.open(filename) as f:
        depth = 0
        for line in f:
            depth += 1
    return depth

snakemake.params.scaling = ''
if hasattr(snakemake.input, 'control'):
    chip_depth = sum([get_depth(x) for x in snakemake.input.chip])
    input_depth = sum([get_depth(x) for x in snakemake.input.control])
    if chip_depth > input_depth:
        snakemake.params.scaling = '--to-large'

#Make sure Macs2 directory exists
subprocess.run(['mkdir', '-p', 'Macs2'])

peaktype = peak_type(snakemake.wildcards.library)
command = ['macs2', 'callpeak', '-t'] +snakemake.input.chip
if hasattr(snakemake.input, 'control'):
    command += ['-c', snakemake.input.control]
#if snakemake.wildcards.library.startswith('DNase') or snakemake.wildcards.library.startswith('ATAC'):
command += ['-f', 'BAM']
#else:
#    command += ['-f', 'BED']
command += ['-n', 'Macs2/{}'.format(snakemake.wildcards.library), '-g', snakemake.config["genomesize"], '--keep-dup', 'all', '--extsize', '200', '-B', '--SPMR']
if snakemake.wildcards.library.startswith('ATAC'):
    command += ['--nomodel', '--shift', '-100']
else:
    command += ['--fix-bimodal']
if snakemake.params.scaling:
    command.append(snakemake.params.scaling)
if peaktype == 'narrow':
    command += ['-q', '0.01']#, '--call-summits']
elif peaktype == 'broad':
    command += ['-q', '0.05', '--broad']
print(command)
subprocess.run([str(arg) for arg in command])

peakfile = 'Macs2/{}_peaks.{}Peak'.format(snakemake.wildcards.library, peaktype)
subprocess.run(['cp', peakfile, snakemake.output.peaks])
