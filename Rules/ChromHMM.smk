
#def model_inputs(wildcards):
    #for assay in PEAK_ASSAYS:
        #if wildcards.type == 'Reads':
            #for library in libraries(assay=assay):
                #if wildcards.type == 'Reads':
                    #yield 'Aligned_Reads/{library}.tagAlign.gz'.format(library=library)
                    #if control_library(library, format='tagAlign'):
                        #yield control_library(library, format='tagAlign')
        #elif wildcards.type == 'Bed':
            #for tissue in tissues_for_assay(assay):
                #yield 'Peak_Calls/{assay}_{tissue}_Combined_Peaks.bed'.format(assay=assay, tissue=tissue)
        #elif wildcards.type == 'IDR':
            #for tissue in tissues_for_assay(assay):
                #yield 'Peak_Calls/{assay}_{tissue}_IDR.{peaktype}Peak'.format(assay=assay, tissue=tissue, peaktype=peak_type(assay))
#
#def individual_model_inputs(wildcards):
    #for assay in PEAK_ASSAYS:
        #for library in libraries(assay=assay,tissue=wildcards.tissue):
            #yield 'Aligned_Reads/{library}.tagAlign.gz'.format(library=library)
            #if control_library(library, format='tagAlign'):
                #yield control_library(library, format='tagAlign')
    

#def training_inputs(wildcards):
    #for data in model_inputs(wildcards):
        #if data not in config['ChromHMM_omit_training']:
            #yield data

def model_inputs(wildcards):
    if wildcards.tissue == 'Joint':
        wildcards.tissue = None
    for library in libraries(skip_input=False, tissue=wildcards.tissue):
        assay, tissue, rep = library.split("_")
        if assay in PEAK_ASSAYS or assay in config['inputs']:# and (assay not in config['no_input'] or wildcards.type=='InputTrack'): 
            yield 'Aligned_Reads/{library}.tagAlign.gz'.format(library=library)
        #elif assay in config['inputs']:
            #yield 'Aligned_Reads/{library}.tagAlign.gz'.format(library=library)

def model_replicate_inputs(wildcards):
    for library in libraries(tissue=wildcards.tissue, rep=wildcards.rep):
        assay, tissue, rep = library.split("_")
        if assay in PEAK_ASSAYS or assay in config['inputs']:# and (assay not in config['no_input'] or wildcards.type=='InputTrack'): 
            yield 'Aligned_Reads/{library}.tagAlign.gz'.format(library=library)
        #elif assay in config['inputs']:
            #yield 'Aligned_Reads/{library}.tagAlign.gz'.format(library=library)

wildcard_constraints:
    states = "[0-9]+",
    sample = "[^_]+_[^_]+"

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
        with open(output.marks, 'w') as f:
            if wildcards.tissue == 'Joint':
                wildcards.tissue = None
            for library in libraries(skip_input=False, tissue=wildcards.tissue):
                assay, tissue, rep = library.split("_")
                if wildcards.type == 'NormInput':
                    if assay in PEAK_ASSAYS:# and assay not in config['no_input']: 
                        f.write("{tissue}\t{assay}\tAligned_Reads/{library}.tagAlign.gz\t{control}\n".format(tissue=tissue, assay=assay, library=library, control=control_library(library, format='tagAlign')))
                elif wildcards.type == 'InputTrack':
                    if assay in PEAK_ASSAYS and assay not in config['inputs']:
                        f.write("{tissue}\t{assay}\tAligned_Reads/{library}.tagAlign.gz\n".format(tissue=tissue, assay=assay, library=library))
                    elif assay in config['inputs']:
                        f.write("{tissue}\tControl\tAligned_Reads/{library}.tagAlign.gz\n".format(tissue=tissue, library=library))

rule replicate_marks:
    output:
        marks = 'ChromHMM/Replicate_Marks_{tissue}_{rep}_{type}.txt',
    run:
        with open(output.marks, 'w') as f:
            for library in libraries(tissue=wildcards.tissue, rep=wildcards.rep):
                assay, tissue, rep = library.split("_")
                if wildcards.type == 'NormInput':
                    if assay in PEAK_ASSAYS:# and assay not in config['no_input']:
                        f.write("{tissue}_{rep}\t{assay}\tAligned_Reads/{library}.tagAlign.gz\t{control}\n".format(tissue=tissue, rep=rep, assay=assay, library=library, control=control_library(library, format='tagAlign')))
                elif wildcards.type == 'InputTrack':
                    if assay not in config['inputs']:
                        f.write("{tissue}_{rep}\t{assay}\tAligned_Reads/{library}.tagAlign.gz\n".format(tissue=tissue, assay=assay, library=library,rep=rep))
                    elif assay in config['inputs']:
                        f.write("{tissue}_{rep}\tControl\tAligned_Reads/{library}.tagAlign.gz\n".format(tissue=tissue, library=library, rep=rep))
                    
rule binarize_data:
    input:
        chroms = 'ChromHMM/Chromosome_Lengths.txt',
        marks = 'ChromHMM/Tissue_Marks_{tissue}_{type}.txt',
        inputs = model_inputs
    output:
        directory('ChromHMM/Binarized_Data_{tissue}_{type}')
    shell:
        'java -mx10000M -jar /home/ckern/ChromHMM/ChromHMM.jar BinarizeBed {input.chroms} . {input.marks} {output}'

rule binarize_replicate:
    input:
        chroms = 'ChromHMM/Chromosome_Lengths.txt',
        marks = 'ChromHMM/Replicate_Marks_{tissue}_{rep}_{type}.txt',
        inputs = model_replicate_inputs
    output:
        directory('ChromHMM/Binarized_Replicate_{tissue}_{rep}_{type}')
    shell:
        'java -mx10000M -jar /home/ckern/ChromHMM/ChromHMM.jar BinarizeBed {input.chroms} . {input.marks} {output}'
        
rule learn_model:
    input:
        chroms = 'ChromHMM/Chromosome_Lengths.txt',
        marks = 'ChromHMM/Tissue_Marks_{tissue}_{type}.txt',
        bindir = 'ChromHMM/Binarized_Data_{tissue}_{type}'
    output:
        emissions = 'ChromHMM/Model_{tissue}_{type}_{states}/emissions_{states}.txt'
    params:
        outdir = 'ChromHMM/Model_{tissue}_{type}_{states}'
    threads: 24
    shell:
        'java -mx10000M -jar /home/ckern/ChromHMM/ChromHMM.jar LearnModel -printposterior -r 300 -p {threads} -l {input.chroms} {input.bindir} {params.outdir} {wildcards.states} {config[ChromHMM_genome]}'
        
rule replicate_segmentation:
    input:
        chroms = 'ChromHMM/Chromosome_Lengths.txt',
        marks = 'ChromHMM/Replicate_Marks_{tissue}_{rep}_{type}.txt',
        bindir = 'ChromHMM/Binarized_Replicate_{tissue}_{rep}_{type}',
        model = 'ChromHMM/Model_Joint_{type}_{states}/model_{states}.txt',
        modeldir = 'ChromHMM/Model_Joint_{type}_{states}'
    output:
        'ChromHMM/Model_Joint_{type}_{states}/{tissue}_{rep}_{states}_segments.bed'
    shell:
        'java -mx10000M -jar /home/ckern/ChromHMM/ChromHMM.jar MakeSegmentation -printposterior {input.model} {input.bindir} {input.modeldir}'
        

##################
# Evaluate Model #
##################

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
        yield 'ChromHMM/Model_{scope}_{type}_{states}/{tissue}_{states}_segments.bed'.format(scope=wildcards.scope, type=wildcards.type, states=wildcards.states, tissue=tissue)

rule get_state_locs_for_clustering:
    input:
        segments = all_chromatin_annotations
    output:
        bed = 'ChromHMM/Model_{scope}_{type}_{states}/State_{state}_Locations.bed'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        'rm -f temp temp2; for file in {input.segments}; do grep {wildcards.state} $file | cut -f1,2,3 >> temp; done; sort -k1,1 -k2,2n temp > temp2;'
        'bedtools merge -i temp2 > {output.bed}; rm -f temp temp2'
        
rule cluster_states:
    input:
        states = rules.get_state_locs_for_clustering.output.bed
    output:
        txt = 'ChromHMM/Model_{scope}_{type}_{states}/State_{state}_Clusters.txt',
        png = 'ChromHMM/Model_{scope}_{type}_{states}/State_{state}_Clusters_Heatmap.png'
    conda:
        '../Envs/python.yaml'
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
