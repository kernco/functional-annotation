
ruleorder: assay_boxplot > assay_boxplot_combined

rule make_library_stats:
    input:
        json = 'Metrics/{library}_Alignment_Stats.json',
        peaks = 'Peak_Calls/{library}_Peaks.bed',
        frip = 'Metrics/{library}_FRiP.txt'
    output:
        json = 'Metrics/{library}_Library_Stats.json'
    script:
        '../Scripts/Get_Library_Stats.py'

rule signal_depth_table:
    input:
        stats = expand('Metrics/{library}_Library_Stats.json', library=libraries(peak_assay_only=True))
    output:
        txt = 'Tables/Signal_Depth.txt',
        csv = 'Tables/Signal_Depth.csv'
    script:
        '../Scripts/Signal_Depth_Table.py'

rule combined_signal_depth_table:
    input:
        stats = ['Metrics/{}_{}_Merged_FRiP.txt'.format(assay, tissue) for assay in PEAK_ASSAYS for tissue in tissues_for_assay(assay)]
    output:
        txt = 'Tables/Merged_Signal_Depth.txt',
        csv = 'Tables/Merged_Signal_Depth.csv'
    script:
        '../Scripts/Merged_Signal_Table.py'

rule make_peak_summary_report:
    input: 
        lambda wildcards: expand(['Peak_Calls/{assay}_{tissue}_Combined_Peaks.bed', 'Metrics/{assay}_{tissue}_FRiP.txt'], assay=wildcards.assay, tissue=LIBRARIES[wildcards.assay]),
        lambda wildcards: expand('Metrics/{library}_FRiP.txt', library=libraries(assay=wildcards.assay))
    output: 
        txt = 'Tables/{assay}_Peak_Summary.txt', 
        csv = 'Tables/{assay}_Peak_Summary.csv'
    params: 
        tissues = lambda wildcards: list(LIBRARIES[wildcards.assay].keys()), 
        reps = lambda wildcards: replicates_for_assay(wildcards.assay)
    script: 
        '../Scripts/Peak_Call_Summary.py'

rule make_idr_peak_summary_report:
    input:
        lambda wildcards: expand(['Peak_Calls/{assay}_{tissue}_IDR.{peaktype}Peak', 'Metrics/{assay}_{tissue}_IDR_FRiP.txt'], assay=wildcards.assay, tissue=[t for t in LIBRARIES[wildcards.assay] if len(list(libraries(assay=wildcards.assay, tissue=t))) > 1], peaktype=peak_type(wildcards.assay))
    output:
        txt = 'Tables/{assay}_IDR_Peak_Summary.txt'
    params:
        tissues = lambda wildcards: list(LIBRARIES[wildcards.assay].keys()), 
        reps = lambda wildcards: replicates_for_assay(wildcards.assay),
        peaktype = lambda wildcards: peak_type(wildcards.assay)
    script:
        '../Scripts/IDR_Summary.py'

######################
# Correlation Graphs #
######################
ruleorder: jaccard_matrix_idr > jaccard_matrix

rule calculate_jaccard:
    input: peaks = lambda wildcards: expand('Peak_Calls/{library}_Peaks.bed', library=libraries(assay=wildcards.assay))
    output: txt='Metrics/{assay}_Jaccard.txt'
    params: hide_assay=True
    conda: '../Envs/deeptools.yaml'
    script: '../Scripts/Make_Jaccard_Matrix.py'

rule jaccard_matrix:
    input: 'Metrics/{assay}_Jaccard.txt'
    output: 'Figures/{assay}_Peak_Similarity.png'
    conda: '../Envs/r.yaml'
    script: '../Scripts/Jaccard_Graph.R'

rule normr_quality_enrichment:
    input:
        treatment = 'Aligned_Reads/{library}.bam',
        control = 'Aligned_Reads/Input_ControlLeghornBursa_1020.bam',#lambda wildcards: control_library("{library}".format(library=wildcards.library)),
        lengths = 'ChromHMM/Chromosome_Lengths.txt'
    output:
        regions = 'Enriched_Regions/{library}_1000fdr0.01.bed'
    params:
        binsize = 1000,
        fdr = 0.01
    conda:
        '../Envs/normr.yaml'
    script:
        '../Scripts/NormR_Enriched_Regions.R'

rule calculate_enriched_regions_jaccard:
    input: peaks=lambda wildcards: expand('Enriched_Regions/{library}_1000fdr0.01.bed', library=libraries(assay=wildcards.assay))
    output: txt='Metrics/{assay}_Enriched_Regions_Jaccard.txt'
    params: hide_assay=True
    conda: '../Envs/deeptools.yaml'
    script: '../Scripts/Make_Jaccard_Matrix.py'

rule enriched_regions_jaccard_graph:
    input: 'Metrics/{assay}_Enriched_Regions_Jaccard.txt'
    output: 'Figures/{assay}_Enriched_Region_Similarity.png'
    conda: '../Envs/r.yaml'
    script: '../Scripts/Jaccard_Graph.R'

rule calculate_all_jaccard:
    input: peaks=lambda wildcards: expand('Peak_Calls/{library}_Peaks.bed', library=peak_libraries())
    output: txt='Metrics/All_Jaccard.txt'
    params: hide_assay=False
    script: '../Scripts/Make_Jaccard_Matrix.py'

rule calculate_jaccard_idr:
    input: peaks=lambda wildcards: expand('IDR/{library}_peaks.{peaktype}Peak', library=libraries(assay=wildcards.prefix.split('_')[0]), peaktype=peak_type(wildcards.prefix.split('_')[0])) + 'IDR/{prefix}.PooledInPr1AndPr2.{peaktype}Peak IDR/{prefix}.PooledInRep1AndRep2.{peaktype}Peak IDR/{prefix}.Pooled_peaks.{peaktype}Peak IDR/{prefix}.Poolr1_peaks.{peaktype}Peak IDR/{prefix}.Poolr2_peaks.{peaktype}Peak Peak_Calls/{prefix}_IDR.{peaktype}Peak'.format(prefix=wildcards.prefix, peaktype=peak_type(wildcards.prefix.split('_')[0])).split()
    output: txt='Metrics/{prefix}_IDR_Jaccard.txt'
    script: '../Scripts/Make_Jaccard_Matrix.py'

rule jaccard_matrix_idr:
    input: 'Metrics/{prefix}_IDR_Jaccard.txt'
    output: 'Figures/{prefix}_IDR_Peak_Similarity.png'
    script: '../Scripts/Jaccard_Graph.R'

rule depth_bedgraph:
    input:
        'Aligned_Reads/{library}.bam'
    output:
        temp('DeepTools/{library}_Depth.bg')
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell:
        'bamCoverage -b {input} -o {output} -of bedgraph -p {threads} -bs 1 -e 200'

rule merge_all_peaks:
    input:
        peaks = lambda wildcards: expand('Peak_Calls/{library}_Peaks.bed', library=libraries(assay=wildcards.assay))
    output:
        'Peak_Calls/{assay}_Merged_Peaks.bed'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        'cat {input.peaks} | sort -k1,1 -k2,2n | bedtools merge -i stdin > {output}'

rule depth_at_peaks:
    input:
        depth = 'DeepTools/{library}_Depth.bg',
        peaks = lambda wildcards: 'Peak_Calls/{assay}_Merged_Peaks.bed'.format(assay=wildcards.library.split('_')[0])
    output:
        'DeepTools/{library}_Peak_Depth.bg'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        'bedtools intersect -a {input.depth} -b {input.peaks} | sort -k1,1 -k2,2n > {output}'
        
rule weighted_jaccard_matrix:
    input: 
        peaks=lambda wildcards: expand('Peak_Calls/{library}_Peaks.bed', library=libraries(assay=wildcards.assay)),
        bedgraphs=lambda wildcards: expand('DeepTools/{library}_Peak_Depth.bg', library=libraries(assay=wildcards.assay))
    output: txt='Metrics/{assay}_Weighted_Jaccard.txt'
    params: hide_assay=True
    conda: '../Envs/bedtools.yaml'
    threads: 24
    script: '../Scripts/Make_Weighted_Jaccard_Matrix.py'

rule weighted_jaccard_graph:
    input: 'Metrics/{assay}_Weighted_Jaccard.txt'
    output: 'Figures/{assay}_Weighted_Jaccard.png'
    conda: '../Envs/r.yaml'
    script: '../Scripts/Jaccard_Graph.R'

################################
# Correlate Peaks with RNA-seq #
################################
rule overlap_genes_with_peaks:
    input: 
        genes = lambda wildcards: 'StringTie/RNASeq_{}_Merged.gtf'.format('_'.join(wildcards.prefix.split('_')[-2:])), 
        peaks = 'Peak_Calls/{prefix}_Peaks.bed'
    output: 
        'Gene_Groups/{prefix}_Pos_Genes.gtf'
    shell: 
        "if [ -s {input.peaks} ]; then bedtools intersect -a {input.genes} -b {input.peaks} -u | grep '\stranscript\s' > {output}; else touch {output}; fi"

rule overlap_genes_without_peaks:
    input: 
        genes = lambda wildcards: 'StringTie/RNASeq_{}_Merged.gtf'.format('_'.join(wildcards.prefix.split('_')[-2:])), 
        peaks = 'Peak_Calls/{prefix}_Peaks.bed'
    output: 
        'Gene_Groups/{prefix}_Neg_Genes.gtf'
    shell: 
        "if [ -s {input.peaks} ]; then bedtools intersect -a {input.genes} -b {input.peaks} -v | grep '\stranscript\s' > {output}; else touch {output}; fi"

rule assay_boxplot:
    input: 
        lambda wildcards: expand(['Gene_Groups/{library}_Pos_Genes.gtf', 'Gene_Groups/{library}_Neg_Genes.gtf'], library=libraries(assay=wildcards.assay))
    output: 
        png='Figures/{assay}_Gene_Expression_Boxplot.png'
    params: 
        labels=lambda wildcards: list(expand(['{library} Pos', '{library} Neg'], library=[' '.join(x.split('_')[1:]) for x in libraries(assay=wildcards.assay)]))
    script: 
        '../Scripts/Expression_Boxplot.py'

rule assay_boxplot_combined:
    input: 
        lambda wildcards: expand(['Gene_Groups/{tissue}_Pos_Genes.gtf', 'Gene_Groups/{tissue}_Neg_Genes.gtf'], tissue=tissues_for_assay(assay=wildcards.prefix.split('_')[0]))
    output: 
        png='Figures/{prefix}_Gene_Expression_Boxplot.png'
    params: 
        labels=lambda wildcards: list(expand(['{tissue} Pos', '{tissue} Neg'], tissue=[' '.join(x.split('_')[1:]) for x in tissues_for_assay(assay=wildcards.prefix.split('_')[0])]))
    script: 
        '../Scripts/Expression_Boxplot.py'
