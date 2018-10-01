from collections import defaultdict

#shell.prefix("module load bio3; ")

####################################
# Load config and detect libraries #
####################################
RAW_FILES, = glob_wildcards("Raw_Reads/{file}.fq.gz")
LIBRARIES = defaultdict(lambda : defaultdict(set))
for filename in RAW_FILES:
    assay, tissue, replicate = filename.split("_")[:3]
    if tissue in config['tissues']:# and replicate in config['reps']:
        LIBRARIES[assay][tissue].add(replicate)
PEAK_ASSAYS = config['broad_peaks'] + config['narrow_peaks']
if 'no_input' not in config:
    config['no_input'] = []
if 'RNASeq' not in config['no_input']:
    config['no_input'].append('RNASeq')
INPUT_ASSAYS = [x for x in PEAK_ASSAYS if x not in config['no_input']]


def libraries(assay=None, tissue=None, rep=None, sample=None, has_input=None, skip_input=True, omit_training=False):
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
        for _tissue in sorted(LIBRARIES[_assay]):
            for _rep in sorted(LIBRARIES[_assay][_tissue]):
                if (not assay or assay == _assay) and (not tissue or tissue == _tissue) and (not rep or rep == _rep):
                    lib = '{assay}_{tissue}_{rep}'.format(assay=_assay, tissue=_tissue, rep=_rep)
                    if omit_training and 'ChromHMM_training_omit' in config and lib in config['ChromHMM_training_omit']:
                        continue
                    yield lib

def peak_libraries():
    for assay in PEAK_ASSAYS:
        for library in libraries(assay=assay):
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
        suffix = 'bam'
    elif format == 'tagAlign':
        suffix = 'tagAlign.gz'
    if assay in config['no_input']:
        return ''
    elif 'override_input' in config and library in config['override_input']:
        return "Aligned_Reads/{}.{}".format(config['override_input'][library], suffix)
    return "Aligned_Reads/{control}_{tissue}_{rep}.{suffix}".format(control=config[assay + "_input"], tissue=tissue, rep=replicate, suffix=suffix)

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

def results(wildcards):
    """ Constructs a list of output files for the pipeline """
    files = []
    files.extend(expand(['Tables/{assay}_Alignment_Summary.txt', 
                     'Figures/{assay}_Pearson_Correlation.png',
                     'Figures/{assay}_Spearman_Correlation.png'],
                     assay=LIBRARIES))
    files.extend(expand(['Tables/{assay}_Quality_Metrics.txt'], assay=INPUT_ASSAYS))
    files.extend(expand(['Tables/{assay}_Peak_Summary.txt', 
                     #'Tables/{assay}_IDR_Peak_Summary.txt',
                     'Figures/{assay}_Peak_Similarity.png'],
                     assay=PEAK_ASSAYS))
    files.extend(expand(['Figures/{assay}_vs_Input_Pearson_Correlation.png', 
                         'Figures/{assay}_vs_Input_Spearman_Correlation.png',
                         'Metrics/{assay}_Coverage.png',
                         'Metrics/{assay}_Fingerprint.png'],
                         assay=[assay for assay in PEAK_ASSAYS if assay not in config['no_input']]))
    files.extend(expand(['Metrics/{sample}_Fingerprint.png', 'Figures/{sample}_TSS_Heatmap.png'],
                         sample=samples()))
    if 'RNASeq' in LIBRARIES:
        files.extend(expand(['Figures/{assay}_Gene_Expression_Boxplot.png'], assay=PEAK_ASSAYS))
    return files

wildcard_constraints:
    assay="[^_]+",
    tissue="[^_]+",
    rep="[^_]+",
    library="[^_]+_[^_]+_[^_]+",
    name="[^\.]+"

rule all:
    input: results

include: "Rules/Align.smk"
include: "Rules/RNAseq.smk"
include: "Rules/ChIPMetrics.smk"
include: "Rules/CallPeaks.smk"
include: "Rules/PeakMetrics.smk"
include: "Rules/ChromHMM.smk"
include: "Rules/ChIPRNACorrelation.smk"
include: "Rules/DeployTrackHub.smk"
include: "Rules/AlleleSpecific.smk"


