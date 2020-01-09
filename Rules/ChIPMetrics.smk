
def deeptools_jsd_inputs(wildcards):
    inputs = {}
    inputs['chip'] = 'Aligned_Reads/{library}.bam'.format(library=wildcards.library),
    inputs['cbai'] = 'Aligned_Reads/{library}.bam.bai'.format(library=wildcards.library)
    if wildcards.library.split('_')[0] not in config['no_input']:
        inputs['control'] = control_library(wildcards.library)
        inputs['cnbai'] = control_library(wildcards.library) + '.bai'
    return inputs

rule fingerprint_sample_graph:
    input:
        bams = lambda wildcards: expand('Aligned_Reads/{library}.bam', library=libraries(sample=wildcards.sample, skip_input=False )),
        bais = lambda wildcards: expand('Aligned_Reads/{library}.bam.bai', library=libraries(sample=wildcards.sample, skip_input=False))
    output:
        'Metrics/{sample}_Fingerprint.png'
    params:
        labels = lambda wildcards: list(libraries(sample=wildcards.sample, skip_input=False))
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell:
        'plotFingerprint -b {input.bams} -plot {output} --labels {params.labels} -p={threads} --extendReads 200 --skipZeros'

rule fingerprint_assay_graph:
    input:
        bams = lambda wildcards: expand('Aligned_Reads/{library}.bam', library=libraries(assay=wildcards.assay, skip_input=False)),
        bais = lambda wildcards: expand('Aligned_Reads/{library}.bam.bai', library=libraries(assay=wildcards.assay, skip_input=False))
    output:
        'Metrics/{assay}_Fingerprint.png'
    params:
        labels = lambda wildcards: list(libraries(assay=wildcards.assay, skip_input=False))
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell:
        'plotFingerprint -b {input.bams} -plot {output} --labels {params.labels} -p={threads} --extendReads 200 --skipZeros'

rule coverage_graph:
    input:
        bams = lambda wildcards: expand('Aligned_Reads/{library}.bam', library=libraries(assay=wildcards.assay, skip_input=False)),
        bais = lambda wildcards: expand('Aligned_Reads/{library}.bam.bai', library=libraries(assay=wildcards.assay, skip_input=False))
    output:
        'Metrics/{assay}_Coverage.png'
    params:
        labels = lambda wildcards: list(libraries(assay=wildcards.assay))
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell:
        'plotCoverage -b {input.bams} -o {output} --labels {params.labels} -p {threads} --extendReads 200'

rule deeptools_jsd:
    input:
        unpack(deeptools_jsd_inputs)
    output:
        metrics = 'Metrics/{library}_DeepTools_Metrics.txt',
        png = temp('Metrics/{library}_Temp.png')
    params: 
        control = lambda wildcards, input: input.control if hasattr(input, 'control') else ''
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell:
        'plotFingerprint -b {input.chip} {params.control} -p={threads} -plot {output.png} --extendReads=200 --skipZeros --outQualityMetrics={output.metrics} --JSDsample={params.control}'

rule bam_coverage:
    input: bam = 'Aligned_Reads/{library}.bam', bai = 'Aligned_Reads/{library}.bam.bai'
    output: 'DeepTools/{library}.bw'
    conda: '../Envs/deeptools.yaml'
    threads: 24
    shell: 'bamCoverage -p={threads} -b {input.bam} -o {output} -of=bigwig --normalizeUsing RPKM --effectiveGenomeSize {config[genomesize]} --ignoreDuplicates --extendReads=200 -bs 20'

rule bam_compare:
    input:
        chipbam = 'Aligned_Reads/{library}.bam',
        chipbai = 'Aligned_Reads/{library}.bam.bai',
        inputbam = lambda wildcards: control_library(wildcards.library),
        inputbai = lambda wildcards: control_library(wildcards.library) + '.bai'
    output: 
        'DeepTools/{library}_vs_Input.bw'
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell: 
        'bamCompare -b1 {input.chipbam} -b2 {input.inputbam} -o {output} -p={threads} --scaleFactorsMethod SES --effectiveGenomeSize {config[genomesize]} --ignoreDuplicates --extendReads=200 -bs 5'

rule multibigwig_summary:
    input: 
        lambda wildcards: expand('DeepTools/{library}.bw', library=libraries(assay=wildcards.assay, skip_input=False))
    output: 
        'DeepTools/{assay}_MultiBigwigSummary.npz'
    params: 
        labels = lambda wildcards: ['"' + ' '.join(x.split('_')[1:]) + '"' for x in list(libraries(assay=wildcards.assay, skip_input=False))]
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell: 
        'multiBigwigSummary bins -p={threads} -b {input} -out {output} --labels {params.labels} --binSize=1000'

rule all_multibigwig_summary:
    input: 
        lambda wildcards: expand('DeepTools/{library}.bw', library=peak_libraries())
    output: 
        'DeepTools/All_MultiBigwigSummary.npz'
    params: 
        labels = lambda wildcards: ['"' + ' '.join(x.split('_')) + '"' for x in peak_libraries()]
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell: 
        'multiBigwigSummary bins -p={threads} -b {input} -out {output} --labels {params.labels} --binSize=1000'

rule all_input_multibigwig_summary:
    input: 
        expand('DeepTools/{library}.bw', library=input_libraries())
    output: 
        'DeepTools/All_Input_MultiBigwigSummary.npz'
    params:
        labels = lambda wildcards: ['"' + ' '.join(x.split('_')) + '"' for x in peak_libraries()]
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell: 
        'multiBigwigSummary bins -p={threads} -b {input} -out {output} --labels {params.labels} --binSize=200'

rule multibigwig_vs_input_summary:
    input: 
        lambda wildcards: expand('DeepTools/{library}_vs_Input_ZScores.bw', library=libraries(assay=wildcards.assay))
    output: 
        'DeepTools/{assay}_vs_Input_MultiBigwigSummary.npz'
    params: 
        labels = lambda wildcards: ['"' + ' '.join(x.split('_')[1:]) + '"' for x in list(libraries(assay=wildcards.assay, skip_input=False))]
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell: 
        'multiBigwigSummary bins -p={threads} -b {input} -out {output} --labels {params.labels} --binSize=1000'# --outRawCounts DeepTools/{wildcards.assay}_vs_Input_BinScores.tab'
        
rule correlation_graph:
    input: 
        'DeepTools/{assaycomp}_MultiBigwigSummary.npz'
    output: 
        'Figures/{assaycomp}_{correlation}_Correlation.png'
    params: 
        lowercase = lambda wildcards: wildcards.correlation.lower()
    conda:
        '../Envs/deeptools.yaml'
    threads: 8
    shell: 
        'plotCorrelation -in {input} --corMethod {params.lowercase} --skipZeros --plotTitle "{wildcards.correlation} Correlation of {wildcards.assaycomp} Read Counts" --whatToPlot heatmap --plotNumbers -o {output} --removeOutliers'

colorlist = ["red", "green", "blue", "orange", "cyan", "purple", "brown", "yellow"]
markerlist = ["s", "o", "v", "p", "*", "P", ">", "<", "^", "d"]
def pca_colors(wildcards):
    if wildcards.assaycomp == "All":
        liblist = list(libraries(peak_assay_only=True, skip_input=True))
        colormap = {x:y for x, y in zip(PEAK_ASSAYS, colorlist)}
        return ' '.join([colormap[x.split('_')[0]] for x in liblist])
    else:
        assay = wildcards.assaycomp.split('_')[0]
        liblist = list(libraries(assay=assay, skip_input=False))
        colormap = {x:y for x, y in zip(tissues_for_assay(assay), colorlist)}
        return ' '.join([colormap[x.split('_')[1]] for x in liblist])

def pca_markers(wildcards):
    if wildcards.assaycomp == "All":
        liblist = list(libraries(peak_assay_only=True, skip_input=False))
        tissues = set([x.split('_')[1] for x in liblist])
        markermap = {x:y for x, y in zip(tissues, markerlist)}
        return ' '.join(['"' + markermap[x.split('_')[1]] + '"' for x in liblist])
    else:
        assay = wildcards.assaycomp.split('_')[0]
        liblist = list(libraries(assay=assay, skip_input=False))
        markermap = {x:y for x, y in zip(replicates_for_assay(assay), markerlist)}
        return ' '.join(['"' + markermap[x.split('_')[2]] + '"' for x in liblist])

rule pca_plot:
    input: 
        'DeepTools/{assaycomp}_MultiBigwigSummary.npz'
    output: 
        'Figures/{assaycomp}_PCA_{pc1}_vs_{pc2}.png'
    params:
        colors = lambda wildcards: pca_colors(wildcards),
        markers = lambda wildcards: pca_markers(wildcards)
    conda:
        '../Envs/deeptools.yaml'
    threads: 8
    shell:
        'plotPCA -in {input} -o {output} --plotTitle "PCA of {wildcards.assaycomp}" --colors {params.colors} --markers {params.markers} --ntop 1000 --transpose --PCs {wildcards.pc1} {wildcards.pc2}'

rule make_quality_metrics_report:
    input:
        stats = lambda wildcards: expand('Metrics/{library}_Alignment_Stats.json', library=libraries(assay=wildcards.assay)), 
        #spp = lambda wildcards: expand('Metrics/{library}.spp_stats.txt', library=libraries(assay=wildcards.assay)),
        deeptools = lambda wildcards: expand('Metrics/{library}_DeepTools_Metrics.txt', library=libraries(assay=wildcards.assay))
    output: 
        txt = 'Tables/{assay}_Quality_Metrics.txt', 
        csv = 'Tables/{assay}_Quality_Metrics.csv'
    params: 
        libraries = lambda wildcards: expand('{library}', library=libraries(assay=wildcards.assay))
    script: '../Scripts/Make_Quality_Metrics.py'
