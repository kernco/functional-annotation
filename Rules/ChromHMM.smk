def model_inputs(wildcards):
    if wildcards.tissue == 'Joint':
        wildcards.tissue = None
    for library in libraries(skip_input=False, tissue=wildcards.tissue):
        assay, tissue, rep = library.split("_")
        if wildcards.type == 'PeakCalls':
            if assay in config['ChromHMM_assays']:
                yield 'Peak_Calls/{assay}_{tissue}_Combined_Peaks.bed'.format(assay=assay, tissue=tissue)
        elif wildcards.type == 'NormR':
            if assay in config['ChromHMM_assays']:
                yield 'Enriched_Regions/{assay}_{tissue}_Combined.bed'.format(assay=assay, tissue=tissue)
        elif assay in config['ChromHMM_assays'] or assay in config['inputs']:
            yield 'Aligned_Reads/{library}.tagAlign.gz'.format(library=library)

def normr_table_inputs(wildcards):
    for library in libraries(assay=wildcards.assay):
        assay, tissue, rep = library.split("_")
        if assay in config['ChromHMM_assays']:
            yield 'Enriched_Regions/{assay}_{tissue}_Combined.bed'.format(assay=assay, tissue=tissue)

def model_replicate_inputs(wildcards):
    for library in libraries(tissue=wildcards.tissue, rep=wildcards.rep):
        assay, tissue, rep = library.split("_")
        if assay in config['ChromHMM_assays'] or assay in config['inputs']:# and (assay not in config['no_input'] or wildcards.type=='InputTrack'): 
            yield 'Aligned_Reads/{library}.tagAlign.gz'.format(library=library)
        #elif assay in config['inputs']:
            #yield 'Aligned_Reads/{library}.tagAlign.gz'.format(library=library)

wildcard_constraints:
    states = "[0-9]+",
    state = "E[0-9]+",
    sample = "[^_]+_[^_]+",
    prefix = "[^_]+_[^_]+"

################
# Create Model #
################
rule chromosome_lengths:
    input:
        genome = config["genome"] 
    output:
        chroms = 'ChromHMM/Chromosome_Lengths.txt'
    script:
        '../Scripts/Get_Chromosome_Lengths.py'

rule tissue_marks:
    output:
        marks = 'ChromHMM/Tissue_Marks_{tissue}_{type}.txt'
    run:
        done = set()
        with open(output.marks, 'w') as f:
            if wildcards.tissue == 'Joint':
                wildcards.tissue = None
            for library in libraries(skip_input=False, tissue=wildcards.tissue):
                if library in config['bad_libraries']:
                    continue
                assay, tissue, rep = library.split("_")
                if wildcards.type == 'NormInput':
                    if assay in config['ChromHMM_assays']:# and assay not in config['no_input']: 
                        f.write("{tissue}\t{assay}\tAligned_Reads/{library}.tagAlign.gz\t{control}\n".format(tissue=tissue, assay=assay, library=library, control=control_library(library, format='tagAlign')))
                elif wildcards.type == 'InputTrack':
                    if assay in config['ChromHMM_assays'] and assay not in config['inputs']:
                        f.write("{tissue}\t{assay}\tAligned_Reads/{library}.tagAlign.gz\n".format(tissue=tissue, assay=assay, library=library))
                    elif assay in config['inputs']:
                        f.write("{tissue}\tControl\tAligned_Reads/{library}.tagAlign.gz\n".format(tissue=tissue, library=library))
                elif wildcards.type == 'PeakCalls':
                    if assay in config['ChromHMM_assays'] and tissue+assay not in done:
                        done.add(tissue+assay)
                        f.write('{tissue}\t{assay}\tPeak_Calls/{assay}_{tissue}_Combined_Peaks.bed\n'.format(tissue=tissue, assay=assay))
                elif wildcards.type == 'NormR':
                    if assay in config['ChromHMM_assays'] and tissue+assay not in done:
                        done.add(tissue+assay)
                        f.write('{tissue}\t{assay}\tEnriched_Regions/{assay}_{tissue}_Combined.bed\n'.format(tissue=tissue, assay=assay))

rule replicate_marks:
    output:
        marks = 'ChromHMM/Replicate_Marks_{tissue}_{rep}_{type}.txt',
    run:
        with open(output.marks, 'w') as f:
            for library in libraries(tissue=wildcards.tissue, rep=wildcards.rep):
                assay, tissue, rep = library.split("_")
                if wildcards.type == 'NormInput':
                    if assay in config['ChromHMM_assays']:# and assay not in config['no_input']:
                        f.write("{tissue}_{rep}\t{assay}\tAligned_Reads/{library}.tagAlign.gz\t{control}\n".format(tissue=tissue, rep=rep, assay=assay, library=library, control=control_library(library, format='tagAlign')))
                elif wildcards.type == 'InputTrack':
                    if assay not in config['inputs']:
                        f.write("{tissue}_{rep}\t{assay}\tAligned_Reads/{library}.tagAlign.gz\n".format(tissue=tissue, assay=assay, library=library,rep=rep))
                    elif assay in config['inputs']:
                        f.write("{tissue}_{rep}\tControl\tAligned_Reads/{library}.tagAlign.gz\n".format(tissue=tissue, library=library, rep=rep))
                    
rule normr_enrichment:
    input:
        treatment = 'Aligned_Reads/{library}.bam',
        control = lambda wildcards: control_library("{library}".format(library=wildcards.library)),
        lengths = 'ChromHMM/Chromosome_Lengths.txt'
    output:
        regions = 'Enriched_Regions/{library}.bed'
    params:
        binsize = 200,
        fdr = 0.05
    conda:
        '../Envs/normr.yaml'
    script:
        '../Scripts/NormR_Enriched_Regions.R'

#rule merge_bams:
    ##input:
        #bams = lambda wildcards: expand('Aligned_Reads/{library}.bam', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue, skip_input=False))
    #output:
        #'Aligned_Reads/{assay}_{tissue}_Merged.bam'
    #conda:
        #'../Envs/samtools.yaml'
    #threads: 8
    #shell:
        #'samtools merge -@{threads} {output} {input} && samtools index {output}'

rule treatment_counts:
    input:
        reads = lambda wildcards: expand('Aligned_Reads/{library}.tagAlign.gz', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue)),
        stats = lambda wildcards: expand('Metrics/{library}_Alignment_Stats.json', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue)),
        chromlens = 'ChromHMM/Chromosome_Lengths.txt'
    output:
        counts = 'Enriched_Regions/{assay}_{tissue}_BinCounts.txt'
    script:
        '../Scripts/Count_Reads_Per_Bin.py'

rule control_counts:
    input:
        reads = lambda wildcards: [control_library(x, format='tagAlign') for x in libraries(assay=wildcards.assay, tissue=wildcards.tissue)],
        stats = lambda wildcards: expand('Metrics/{library}_Alignment_Stats.json', library=[control_library(x, format='basename') for x in libraries(assay=wildcards.assay, tissue=wildcards.tissue)]),
        chromlens = 'ChromHMM/Chromosome_Lengths.txt'
    output:
        counts = 'Enriched_Regions/{assay}_{tissue}_Control_BinCounts.txt'
    script:
        '../Scripts/Count_Reads_Per_Bin.py'

rule merged_regimer:
    input:
        treatment = 'Enriched_Regions/{assay}_{tissue}_BinCounts.txt',
        control = 'Enriched_Regions/{assay}_{tissue}_Control_BinCounts.txt',
        lengths = 'ChromHMM/Chromosome_Lengths.txt',
    output:
        regions = 'Enriched_Regions/{assay}_{tissue}_Regimes.bed'
    params:
        binsize = 200,
        fdr = 0.05
    #conda:
        #'../Envs/normr.yaml'
    script:
        '../Scripts/RegimeR_Merged.R'

rule regimer_enrichment:
    input:
        treatment = 'Aligned_Reads/{library}.bam',
        control = lambda wildcards: control_library("{library}".format(library=wildcards.library)),
        lengths = 'ChromHMM/Chromosome_Lengths.txt',
        stats = 'Metrics/{library}_Alignment_Stats.json'
    output:
        regions = 'Enriched_Regions/{library}_Regimes.bed'
    params:
        binsize = 200,
        fdr = 0.05
    conda:
        '../Envs/normr.yaml'
    script:
        '../Scripts/RegimeR_Enriched_Regions.R'

rule split_regimer:
    input:
        infile = 'Enriched_Regions/{root}_Regimes.bed'
    output:
        outfile = 'Enriched_Regions/{root}_Regime{regime}.bed',
    script:
        '../Scripts/Split_Regime.py'

rule combine_regimer:
    input:
        regime1 = lambda wildcards: expand('Enriched_Regions/{library}_Regime1.bed', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue)),
        regime2 = lambda wildcards: expand('Enriched_Regions/{library}_Regime2.bed', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue)),
        merged1 = 'Enriched_Regions/{assay}_{tissue}_Merged_Regime1.bed',
        merged2 = 'Enriched_Regions/{assay}_{tissue}_Merged_Regime2.bed'
    output:
        combinereg1 = "Enriched_Regions/{assay}_{tissue}_Combined_Regime1.bed",
        combinereg2 = "Enriched_Regions/{assay}_{tissue}_Combined_Regime2.bed"
    script:
        '../Scripts/Combine_Regime_Replicates.py'

rule combine_normr:
    input:
        lambda wildcards: expand('Enriched_Regions/{library}.bed', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue))
    output:
        regions = 'Enriched_Regions/{assay}_{tissue}_Combined.bed'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        'bedtools intersect -a {input[0]} -b {input[1]} -u | bedtools merge -i stdin > {output}'

rule normr_summary_table:
    input:
        inputs = normr_table_inputs
    output:
        'Tables/{assay}_NormR_Summary.csv'
    params: 
        tissues = lambda wildcards: list(LIBRARIES[wildcards.assay].keys()), 
        reps = lambda wildcards: replicates_for_assay(wildcards.assay)
    script:
        '../Scripts/NormR_Summary_Table.py'

rule binarize_data:
    input:
        chroms = 'ChromHMM/Chromosome_Lengths.txt',
        marks = 'ChromHMM/Tissue_Marks_{tissue}_{type}.txt',
        inputs = model_inputs
    output:
        directory('ChromHMM/Binarized_Data_{tissue}_{type}')
    conda:
        '../Envs/chromHMM.yaml'
    params:
        peaks = lambda wildcards: '-peaks' if wildcards.type in ['PeakCalls', 'NormR'] else ''
    threads: 12
    shell:
        'ChromHMM.sh BinarizeBed {params.peaks} {input.chroms} . {input.marks} {output}'

rule binarize_replicate:
    input:
        chroms = 'ChromHMM/Chromosome_Lengths.txt',
        marks = 'ChromHMM/Replicate_Marks_{tissue}_{rep}_{type}.txt',
        inputs = model_replicate_inputs
    output:
        'ChromHMM/Binarized_Replicate_{tissue}_{rep}_{type}'
    conda:
        '../Envs/chromHMM.yaml'
    threads: 12
    shell:
        'ChromHMM.sh BinarizeBed {input.chroms} . {input.marks} {output}'
        
rule learn_model:
    input:
        chroms = 'ChromHMM/Chromosome_Lengths.txt',
        marks = 'ChromHMM/Tissue_Marks_{tissue}_{type}.txt',
        bindir = 'ChromHMM/Binarized_Data_{tissue}_{type}'
    output:
        emissions = 'ChromHMM/Model_{tissue}_{type}_{states}/emissions_{states}.txt'
    params:
        outdir = 'ChromHMM/Model_{tissue}_{type}_{states}'
    conda:
        '../Envs/chromHMM.yaml'
    threads: 12
    shell:
        'ChromHMM.sh LearnModel -printposterior -p {threads} -l {input.chroms} {input.bindir} {params.outdir} {wildcards.states} {config[ChromHMM_genome]}'
        
rule replicate_segmentation:
    input:
        chroms = 'ChromHMM/Chromosome_Lengths.txt',
        marks = 'ChromHMM/Replicate_Marks_{tissue}_{rep}_{type}.txt',
        bindir = 'ChromHMM/Binarized_Replicate_{tissue}_{rep}_{type}',
        model = 'ChromHMM/Model_Joint_{type}_{states}/model_{states}.txt',
        modeldir = 'ChromHMM/Model_Joint_{type}_{states}'
    output:
        'ChromHMM/Model_Joint_{type}_{states}/{tissue}_{rep}_{states}_segments.bed'
    conda:
        '../Envs/chromHMM.yaml'
    threads: 12
    shell:
        'ChromHMM.sh MakeSegmentation -printposterior {input.model} {input.bindir} {input.modeldir}'
        

rule split_states:
    input:
        'ChromHMM/{model}/{prefix}_segments.bed'
    output:
        'ChromHMM/{model}/{prefix}_{state}.bed'
    shell:
        'grep {wildcards.state}$ {input} > {output}'

def tissue_specific_inputs(wildcards):
    maintissue, states = wildcards.prefix.split('_')
    inputs = []
    for tissue in config['ChromHMM_tissues']:
        if tissue != maintissue:
            inputs.append('{model}/{tissue}_{states}_{state}.bed'.format(model=wildcards.model, tissue=tissue, states=states, state=wildcards.state))
            if hasattr(wildcards, 'more'):
                more = wildcards.more.split('E')[1:]
                for s in more:
                    inputs.append('{model}/{tissue}_{states}_{state}.bed'.format(model=wildcards.model, tissue=tissue, states=states, state='E' + s))
    return inputs
    
rule tissue_specific_state:
    input:
        specific = '{model}/{prefix}_{state}.bed',
        others = tissue_specific_inputs
    output:
        '{model}/{prefix}_Specific_{state}.bed'
    shell:
        'bedtools intersect -a {input.specific} -b {input.others} -v > {output}'

rule tissue_specific_state_alternate:
    input:
        specific = '{model}/{prefix}_{state}.bed',
        others = tissue_specific_inputs
    output:
        '{model}/{prefix}_Specific_{state}-{more}.bed'
    shell:
        'bedtools intersect -a {input.specific} -b {input.others} -v > {output}'
        
rule get_segment_seqs:
    input:
        bed = '{model}/{spec}_{suffix}.bed',
        genome = lambda wildcards: genomes['{}'.format(wildcards.spec)]
    output:
        '{model}/{spec}_{suffix}.fa'
    shell:
        'bedtools getfasta -fi {input.genome} -bed {input.bed} > {output}'

##################
# Evaluate Model #
##################

rule create_tissue_models:
    input:
        models = lambda wildcards: expand('ChromHMM/Model_{tissue}_{type}_{states}/emissions_{states}.txt', tissue=config['tissues'], type=wildcards.type, states=wildcards.states)
    output:
        modeldir = 'ChromHMM/Tissue_Models_{type}_{states}'
    run:
        shell('mkdir -p {}'.format(output.modeldir))
        for model in input.models:
            tissue = model.split('_')[1]
            shell('cp {} {}/emissions_{tissue}.txt'.format(model, output.modeldir, tissue=tissue))
        
rule test_num_states:
    input:
        testmodel = 'ChromHMM/Model_Joint_{type}_{states}/emissions_{states}.txt',
        tissuemodels = 'ChromHMM/Tissue_Models_{type}_{states}'
    output:
        txt = 'ChromHMM/Correlation_Tests/{type}_{states}_Comparison.txt'
    params:
        prefix = 'ChromHMM/Correlation_Tests/{type}_{states}_Comparison'
    conda:
        '../Envs/chromHMM.yaml'
    shell:
        'ChromHMM.sh CompareModels {input.testmodel} {input.tissuemodels} {params.prefix}'

rule pairwise_overlap:
    input:
        bed1 = 'ChromHMM/Model_{scope}_{type}_{states}/{tissue1}_{states}_segments.bed',
        bed2 = 'ChromHMM/Model_{scope}_{type}_{states}/{tissue2}_{states}_segments.bed'
    output:
        'ChromHMM/Model_{scope}_{type}_{states}/{tissue1}_{tissue2}_Intersect.bed'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        'bedtools intersect -a {input.bed1} -b {input.bed2} -wo > {output}'

rule pairwise_replicate_overlap:
    input:
        bed1 = 'ChromHMM/Model_{scope}_{type}_{states}/{tissue1}_{rep1}_{states}_segments.bed',
        bed2 = 'ChromHMM/Model_{scope}_{type}_{states}/{tissue2}_{rep2}_{states}_segments.bed'
    output:
        'ChromHMM/Model_{scope}_{type}_{states}/{tissue1}_{rep1}_{tissue2}_{rep2}_Intersect.bed'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        'bedtools intersect -a {input.bed1} -b {input.bed2} -wo > {output}'

def pairwise_enrichment_inputs(wildcards):
    tissues = sorted(config['tissues'])
    for i, tissue1 in enumerate(tissues):
        for tissue2 in tissues[i+1:]:
            yield 'ChromHMM/Model_{scope}_{type}_{states}/{tissue1}_{tissue2}_Intersect.bed'.format(scope=wildcards.scope, type=wildcards.type, states=wildcards.states, tissue1=tissue1, tissue2=tissue2)

rule pairwise_enrichment_table:
    input:
        pairwise_enrichment_inputs
    output:
        csv = 'ChromHMM/Model_{scope}_{type}_{states}/Pairwise_Enrichment.csv'
    script:
        '../Scripts/Chromatin_State_Similarity.py'

def across_replicate_enrichment_inputs(wildcards):
    tissues = sorted(config['tissues'])
    for i, tissue1 in enumerate(tissues):
        for tissue2 in tissues[i+1:]:
            for rep in replicates_for_assay(config['narrow_peaks'][0]):
                yield 'ChromHMM/Model_{scope}_{type}_{states}/{tissue1}_{rep}_{tissue2}_{rep}_Intersect.bed'.format(scope=wildcards.scope, type=wildcards.type, states=wildcards.states, tissue1=tissue1, tissue2=tissue2, rep=rep)

def between_replicate_enrichment_inputs(wildcards):
    reps = replicates_for_assay(config['narrow_peaks'][0])
    for i, rep1 in enumerate(reps):
        for rep2 in reps[i+1:]:
            for tissue in config['tissues']:
                yield 'ChromHMM/Model_{scope}_{type}_{states}/{tissue}_{rep1}_{tissue}_{rep2}_Intersect.bed'.format(scope=wildcards.scope, type=wildcards.type, states=wildcards.states, tissue=tissue, rep1=rep1, rep2=rep2)

rule across_replicate_enrichment_table:
    input:
        across_replicate_enrichment_inputs
    output:
        csv = 'ChromHMM/Model_{scope}_{type}_{states}/Across_Replicate_Enrichment.csv'
    script:
        '../Scripts/Chromatin_State_Similarity.py'

rule between_replicate_enrichment_table:
    input:
        between_replicate_enrichment_inputs
    output:
        csv = 'ChromHMM/Model_{scope}_{type}_{states}/Between_Replicate_Enrichment.csv'
    script:
        '../Scripts/Chromatin_State_Similarity.py'
    
def all_chromatin_annotations(wildcards):
    for tissue in config['tissues']:
        yield 'ChromHMM/Model_{scope}_{type}_{states}/{tissue}_{states}_200bp_segments.bed'.format(scope=wildcards.scope, type=wildcards.type, states=wildcards.states, tissue=tissue)

rule chunk_chromatin_annotation:
    input:
        infile = 'ChromHMM/Model_{scope}_{type}_{states}/{tissue}_{states}_segments.bed'
    output:
        outfile = 'ChromHMM/Model_{scope}_{type}_{states}/{tissue}_{states}_200bp_segments.bed'
    script:
        '../Scripts/Chunk_Bed.py'

rule get_state_locs_for_clustering:
    input:
        segments = all_chromatin_annotations
    output:
        bed = 'ChromHMM/Model_{scope}_{type}_{states}/State_{state}_Locations.bed'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        """rm -f temp temp2; for file in {input.segments}; do grep {wildcards.state} $file | awk -v OFS="\\t" '{{print $1, $2+1, $3}}' >> temp; done; sort -k1,1 -k2,2n temp > temp2; """
        'bedtools merge -i temp2 > {output.bed}; rm -f temp temp2'
        
rule cluster_states:
    input:
        states = rules.get_state_locs_for_clustering.output.bed
    output:
        txt = 'ChromHMM/Model_{scope}_{type}_{states}/State_{state}_Clusters.txt',
        png = 'ChromHMM/Model_{scope}_{type}_{states}/State_{state}_Clusters_Heatmap.png'
    conda:
        '../Envs/python.yaml'
    threads: 48
    script:
        '../Scripts/Cluster_Chromatin_States.py'

rule create_tss_gtf:
    input:
        annotation = config['annotation']
    output:
        gtf = 'ChromHMM/TSS.gtf'
    script:
        "../Scripts/Make_TSS_GTF.py"

rule sort_tss_gtf:
    input:
        rules.create_tss_gtf.output.gtf
    output:
        'ChromHMM/TSS_Sorted.gtf'
    shell:
        'sort -k1,1 -k2,2n {input} > {output}'

rule assign_states_to_gene:
    input:
        clusters = rules.cluster_states.output.txt,
        tss = rules.sort_tss_gtf.output
    output:
        txt = 'ChromHMM/Model_{scope}_{type}_{states}/State_{state}_Promoter_Assignment.txt'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        'bedtools closest -a {input.clusters} -b {input.tss} -d > {output}'

rule go_enrichment_analysis:
    input:
        clusters = rules.assign_states_to_gene.output.txt
    output:
        txt = 'ChromHMM/Model_{scope}_{type}_{states}/State_{state}_{state_type}_GO_Enrichment.txt'
    conda:
        '../Envs/python.yaml'
    script:
        '../Scripts/State_Cluster_GO_Enrichment.py'

rule split_cluster_files:
    input:
        clusters = rules.cluster_states.output.txt
    output:
        dir = 'ChromHMM/Model_{scope}_{type}_{states}/State_{state}_Clusters'
    run:
        shell('mkdir -p {output}')
        outfiles = {}
        with open(input.clusters) as f:
            for line in f:
                cols = line.split()
                name = chr(int(cols[-1]) + 65)
                if cols[-1] not in outfiles:
                    outfiles[cols[-1]] = open('ChromHMM/Model_{scope}_{type}_{states}/State_{state}_Clusters/Cluster_{cluster}.bed'.format(states=wildcards.states, state=wildcards.state, cluster=name), 'w')
                outfiles[cols[-1]].write("{}\t{}\t{}\n".format(cols[0], cols[1], int(cols[1]) + 200))

rule state_cluster_motifs:
    input:
        background = rules.get_state_locs_for_clustering.output.bed,
        clusters = 'ChromHMM/Model_{scope}_{type}_{states}/State_{state}_Clusters'
    output:
        homer = 'ChromHMM/Model_{scope}_{type}_{states}/State_{state}_Clusters/Cluster_{cluster}_Motifs'
    conda:
        '../Envs/homer.yaml'
    threads: 12
    shell:
        'findMotifsGenome.pl {input.clusters}/Cluster_{wildcards.cluster}.bed {config[genome]} {output} -size given -bg {input.background} -noweight -p 12'
        
rule assign_states_to_tss:
    input:
        tss = 'StringTie/Merged_TSS.bed',
        model = 'ChromHMM/Model_{scope}_{type}_{states}/{tissue}_{states}_dense.bed'
    output:
        overlap = 'ChromHMM/Model_{scope}_{type}_{states}/{tissue}_TSS_states.txt'
    params:
        segments = lambda wildcards: 'ChromHMM/Model_{scope}_{type}_{states}/{tissue}_{states}_dense.bed'.format(type=wildcards.type, scope=wildcards.scope, states=wildcards.states, tissue=wildcards.tissue)
    shell:
        'bedtools intersect -a {input.tss} -b {params.segments} -wa -wb > {output.overlap}'

rule tss_state_boxplot:
    input:
        overlap = 'ChromHMM/Model_{scope}_{type}_{states}/{tissue}_TSS_states.txt',
        expression = lambda wildcards: ['StringTie/{library}_Expression.tab'.format(library=library) for library in libraries(assay='RNASeq', tissue=wildcards.tissue)]
    output:
        png = 'Figures/ChromHMM_{states}State_{scope}_{type}_{tissue}_TSS_Boxplot.png'
        #roc = 'Figures/ChromHMM_{states}State{type}_{tissue}_ROC.png'
    script:
        '../Scripts/ChromHMM_TSS_Boxplot.py'
