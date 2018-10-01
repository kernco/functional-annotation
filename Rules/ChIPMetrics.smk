
def deeptools_jsd_inputs(wildcards):
    inputs = {}
    inputs['chip'] = 'Aligned_Reads/{library}.bam'.format(library=wildcards.library),
    inputs['cbai'] = 'Aligned_Reads/{library}.bam.bai'.format(library=wildcards.library)
    if wildcards.library.split('_')[0] not in config['no_input']:
        inputs['control'] = control_library(wildcards.library)
        inputs['cnbai'] = control_library(wildcards.library) + '.bai'
    return inputs

rule spp_stats:
    input: 
        'Aligned_Reads/{library}.bam'
    output: 
        stats = 'Metrics/{library}_spp_stats.txt', 
        figure = 'Figures/{library}_Cross_Correlation.pdf'
    threads: 12
    conda:
        '../Envs/r.yaml'
    shell: 
        'Rscript /home/ckern/phantompeakqualtools/run_spp.R -c={input} -rf -out={output.stats} -p={threads} -s=0:2:400 -savp={output.figure}'

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
    output: 'DeepTools/{library}_vs_Input.bw'
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell: 'bamCompare -b1 {input.chipbam} -b2 {input.inputbam} -o {output} -p={threads} --normalizeUsing RPGC --effectiveGenomeSize {config[genomesize]} --ignoreDuplicates --extendReads=200'

rule multibigwig_summary:
    input: lambda wildcards: expand('DeepTools/{library}.bw', library=libraries(assay=wildcards.assay, skip_input=False))
    output: 'DeepTools/{assay}_MultiBigwigSummary.npz'
    params: labels = lambda wildcards: list(libraries(assay=wildcards.assay, skip_input=False))
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell: 'multiBigwigSummary bins -p={threads} -b {input} -out {output} --labels {params.labels} --binSize=1000'

rule multibigwig_vs_input_summary:
    input: lambda wildcards: expand('DeepTools/{library}_vs_Input.bw', library=libraries(assay=wildcards.assay))
    output: 'DeepTools/{assay}_vs_Input_MultiBigwigSummary.npz'
    params: labels = lambda wildcards: list(libraries(assay=wildcards.assay))
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell: 'multiBigwigSummary bins -p={threads} -b {input} -out {output} --labels {params.labels} --binSize=1000'
        
rule correlation_graph:
    input: 'DeepTools/{assaycomp}_MultiBigwigSummary.npz'
    output: 'Figures/{assaycomp}_{correlation}_Correlation.png'
    params: lowercase = lambda wildcards: wildcards.correlation.lower()
    conda:
        '../Envs/deeptools.yaml'
    shell: 'plotCorrelation -in {input} --corMethod {params.lowercase} --skipZeros --plotTitle "{wildcards.correlation} Correlation of {wildcards.assaycomp} Read Counts" --whatToPlot heatmap --plotNumbers -o {output} --removeOutliers'

rule make_quality_metrics_report:
    input:
        bams = lambda wildcards: expand('Aligned_Reads/{library}.bam', library=libraries(assay=wildcards.assay)), 
        spp = lambda wildcards: expand('Metrics/{library}_spp_stats.txt', library=libraries(assay=wildcards.assay)),
        deeptools = lambda wildcards: expand('Metrics/{library}_DeepTools_Metrics.txt', library=libraries(assay=wildcards.assay))
    output: 
        txt = 'Tables/{assay}_Quality_Metrics.txt', 
        csv = 'Tables/{assay}_Quality_Metrics.csv'
    params: 
        libraries = lambda wildcards: expand('{library}', library=libraries(assay=wildcards.assay))
    script: '../Scripts/Make_Quality_Metrics.py'
