from collections import defaultdict

#shell.prefix("module load bio3; ")

####################################
# Load config and detect libraries #
####################################
RAW_FILES, = glob_wildcards("Raw_Reads/{file}.fq.gz")
LIBRARIES = defaultdict(lambda : defaultdict(set))
PAIRED_LIBS = set()
ASSAYS = set()
for filename in RAW_FILES:
    assay, tissue, replicate = filename.split("_")[:3]
    if assay == 'RRBS':
        continue
    ASSAYS.add(assay)
    if filename.endswith("_R1"):
        PAIRED_LIBS.add("{}_{}_{}".format(assay, tissue, replicate))
    if tissue in config['tissues']:# and replicate in config['reps']:
        LIBRARIES[assay][tissue].add(replicate)
PEAK_ASSAYS = config['broad_peaks'] + config['narrow_peaks']
if 'no_input' not in config:
    config['no_input'] = []
if 'RNASeq' not in config['no_input']:
    config['no_input'].append('RNASeq')
INPUT_ASSAYS = [x for x in PEAK_ASSAYS if x not in config['no_input']]


def libraries(assay=None, tissue=None, rep=None, sample=None, has_input=None, skip_input=True, omit_training=False, peak_assay_only=False):
    """ Returns a list of libraries matching the given parameters. """
    if sample:
        tissue, rep = sample.split("_")
    for _assay in sorted(LIBRARIES):
        if skip_input and _assay in config['inputs']:
            continue
        if has_input and _assay in config['no_input']:
            continue
        if has_input == False and _assay not in config['no_input']:
            continue
        if peak_assay_only and _assay not in PEAK_ASSAYS:
            continue
        for _tissue in sorted(LIBRARIES[_assay]):
            for _rep in sorted(LIBRARIES[_assay][_tissue]):
                if (not assay or assay == _assay) and (not tissue or tissue == _tissue) and (not rep or rep == _rep):
                    lib = '{assay}_{tissue}_{rep}'.format(assay=_assay, tissue=_tissue, rep=_rep)
                    if omit_training and 'ChromHMM_training_omit' in config and lib in config['ChromHMM_training_omit']:
                        continue
                    yield lib

def library(**kwargs):
    return list(libraries(**kwargs))[0]

def peak_libraries():
    for assay in PEAK_ASSAYS:
        for library in libraries(assay=assay):
            yield library

def input_libraries():
    for library in libraries(skip_input=False):
        assay, group, rep = library.split('_')
        if assay in config['inputs']:
            yield library

def tissues_for_assay(assay):
    return sorted(LIBRARIES[assay])

def samples(assay=None):
    if assay:
        for tissue in LIBRARIES[assay]:
            for rep in LIBRARIES[assay][tissue]:
                yield "{}_{}".format(tissue, rep)
    else:
        for assay in LIBRARIES:
            for sample in samples(assay):
                yield sample

def combined_libraries(assay=None):
    if isinstance(assay, list):
        for each in assay:
            yield from combined_libraries(each)
    elif assay:
        for tissue in LIBRARIES[assay]:
            yield '{}_{}'.format(assay, tissue)
    else:
        for assay in LIBRARIES:
            yield from combined_libraries(assay)

def replicates_for_assay(assay):
    """ Given an assay, returns a list of all replicates with libraries of that assay """
    reps = set()
    for tissue in LIBRARIES[assay]:
        reps.update(LIBRARIES[assay][tissue])
    return sorted(list(reps))

def other_replicate(library):
    """ Given a library, returns the library of the same assay/tissue
        from the other replicate """
    assay, tissue, replicate = library.split("_")
    for rep in sorted(LIBRARIES[assay][tissue]):
        if rep != replicate:
            return "{assay}_{tissue}_{rep}".format(assay=assay, tissue=tissue, rep=rep)
        return library #failsafe for 1 replicate

def control_library(library, format='BAM'):
    """ Given a ChIP library, returns the corresponding input/control library """
    assay, tissue, replicate = library.split("_")
    if format == 'BAM':
        prefix = 'Aligned_Reads/'
        suffix = '.bam'
    elif format == 'tagAlign' or format == 'tagAlign.gz':
        prefix = 'Aligned_Reads/'
        suffix = '.tagAlign.gz'
    elif format == 'basename':
        prefix = ''
        suffix = ''
    if assay in config['no_input']:
        return ''
    elif 'override_input' in config and library in config['override_input']:
        return prefix + config['override_input'][library] + suffix
    elif library.split('_')[-1] == 'Merged':
        return prefix + library + 'Input' + suffix
    return prefix + '{control}_{tissue}_{rep}'.format(control=config[assay + "_input"], tissue=tissue, rep=replicate) + suffix

def peak_type(library):
    """ Given a library or assay, returns whether the assay is 'broad' or 'narrow' """
    if "_" in library:
        assay, tissue, replicate = library.split("_")
    else:
        assay = library
    if assay in config['broad_peaks']:
        return 'broad'
    elif assay in config['narrow_peaks']:
        return 'narrow'
    else:
        raise Exception

def tables(wildcards):
    files = []
    files.extend(expand('Tables/{assay}_Alignment_Summary.txt', assay=[x for x in ASSAYS if x != "RRBS"]))
    files.extend(expand('Tables/{assay}_Quality_Metrics.txt', assay=INPUT_ASSAYS))
    files.extend(expand('Tables/{assay}_Peak_Summary.txt', assay=PEAK_ASSAYS))
    files.append('Tables/Signal_Depth.txt')
    files.append('Tables/Merged_Signal_Depth.txt')
    if 'RNASeq' in LIBRARIES:
        files.append('Tables/RNASeq_Alignment_Summary.txt')
        #files.append('Tables/Gene_Expression_EdgeR_Report.csv')
    return files

def figures(wildcards):
    """ Constructs a list of output files for the pipeline """
    files = []
    files.extend(expand(['Figures/{assay}_Pearson_Correlation.png',
                     'Figures/{assay}_Spearman_Correlation.png', 
                     'Figures/{assay}_PCA_1_vs_2.png'],
                     assay=PEAK_ASSAYS))
    files.extend(expand('Figures/{assay}_Peak_Similarity.png', assay=PEAK_ASSAYS))
    files.extend(expand(['Figures/{assay}_vs_Input_Pearson_Correlation.png', 
                         'Figures/{assay}_vs_Input_Spearman_Correlation.png',
                         'Figures/{assay}_vs_Input_PCA_1_vs_2.png',
                         'Figures/{assay}_vs_Input_PCA_1_vs_3.png',
                         'Figures/{assay}_vs_Input_PCA_2_vs_3.png',
                         'Metrics/{assay}_Coverage.png',
                         'Metrics/{assay}_Fingerprint.png'],
                         assay=[assay for assay in PEAK_ASSAYS if assay not in config['no_input']]))
    files.extend(expand(['Metrics/{sample}_Fingerprint.png'],
                         sample=samples()))
    if 'RNASeq' in LIBRARIES:
        #files.extend(expand(['Figures/{sample}_TSS_Heatmap.png'], sample=samples(assay='RNASeq')))
        files.append('Figures/TPM_Density_Plot.png')
        files.append('Figures/TSI_Density_Plot.png')
    return files

def results(wildcards):
    return tables(wildcards) + figures(wildcards) + ['Track_Hub/trackDb.txt']

wildcard_constraints:
    assay="[^_/]+",
    tissue="[^_/]+",
    rep="[^_/]+",
    library="[^_/]+_[^_/]+_[^_/]+"
    #name="[^\.]+"

rule all:
    input: results

rule Make_Tables:
    input: tables

rule input_table:
    output: "Tables/Input_Library_Mapping.txt"
    run:
        outfile = open(output[0], 'w')
        for library in libraries(skip_input=True):
            inputlib = control_library(library, format='basename')
            if inputlib:
                print("{}\t{}".format(library, inputlib), file=outfile)

include: "Rules/Align.smk"
include: "Rules/RNAseq.smk"
include: "Rules/ChIPMetrics.smk"
include: "Rules/CallPeaks.smk"
include: "Rules/IDR.smk"
include: "Rules/PeakMetrics.smk"
include: "Rules/ChromHMM.smk"
include: "Rules/ChIPRNACorrelation.smk"
include: "Rules/DeployTrackHub.smk"
include: "Rules/AlleleSpecific.smk"
include: "Rules/Methylation.smk"


