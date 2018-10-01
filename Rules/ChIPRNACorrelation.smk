def matrix_bigwigs(wildcards):
    bigwigs = expand('DeepTools/{library}_vs_Input.bw', library=libraries(assay='H3K4me3', tissue=wildcards.tissue, rep=wildcards.rep))
    bigwigs += expand('DeepTools/{library}_vs_Input.bw', library=libraries(assay='H3K27me3', tissue=wildcards.tissue, rep=wildcards.rep))
    bigwigs += expand('DeepTools/{library}_vs_Input.bw', library=libraries(assay='H3K4me1', tissue=wildcards.tissue, rep=wildcards.rep))
    bigwigs += expand('DeepTools/{library}_vs_Input.bw', library=libraries(assay='H3K27ac', tissue=wildcards.tissue, rep=wildcards.rep))
    bigwigs += expand('DeepTools/{library}.bw', library=libraries(assay='DNaseSeq', tissue=wildcards.tissue, rep=wildcards.rep))
    #bigwigs += expand('DeepTools/{library}.bw', library=libraries(tissue=wildcards.tissue, rep=wildcards.rep, has_input=False))
    return bigwigs

ruleorder: rna_bam_coverage > bam_coverage

rule rna_bam_coverage:
    input:
        bam = 'Aligned_Reads/RNASeq_{sample}.bam',
        bai = 'Aligned_Reads/RNASeq_{sample}.bam.bai'
    output:
        'DeepTools/RNASeq_{sample}.bw'
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell:
        'bamCoverage -p {threads} -b {input.bam} -o {output} -of=bigwig --normalizeUsing RPGC --effectiveGenomeSize {config[genomesize]}'

rule compute_matrix:
    input:
        bigwigs = matrix_bigwigs,
        genes = config['annotation']
    output:
        'DeepTools/{tissue}_{rep}_matrix.mat.gz'
    params:
        labels = lambda wildcards: [x.split("/")[1].split("_")[0] for x in matrix_bigwigs(wildcards)]
    threads: 24
    conda:
        '../Envs/deeptools.yaml'
    shell:
        'computeMatrix scale-regions -S {input.bigwigs} -R {input.genes} --beforeRegionStartLength 2000 --afterRegionStartLength 2000 -out {output} -m 2000 -p {threads} --transcriptID transcript --samplesLabel {params.labels}'

rule compute_matrix_ref:
    input:
        bigwigs = matrix_bigwigs,
        genes = config['annotation']
    output:
        'DeepTools/{tissue}_{rep}_matrix_ref.mat.gz'
    params:
        labels = lambda wildcards: [x.split("/")[1].split("_")[0] for x in matrix_bigwigs(wildcards)]
    threads: 24
    conda:
        '../Envs/deeptools.yaml'
    shell:
        'computeMatrix reference-point -S {input.bigwigs} -R {input.genes} --beforeRegionStartLength 3000 --afterRegionStartLength 3000 -out {output} -p {threads} --transcriptID transcript --samplesLabel {params.labels}'

rule make_tss_region_bed:
    input:
        tab = 'StringTie/RNASeq_{tissue}_{rep}_Annotation_Expression.tab',
        chromsizes = config['chromsizes']
    output:
        bed = 'DeepTools/{tissue}_{rep}_Sorted_TSS.bed'
    script:
        '../Scripts/Gene_TSS_Region_Bed.py'

rule compute_matrix_annotation:
    input:
        bigwigs = matrix_bigwigs,
        genes = rules.make_tss_region_bed.output
    output:
        'DeepTools/{tissue}_{rep}_matrix_annotation.mat.gz'
    params:
        labels = lambda wildcards: [x.split("/")[1].split("_")[0] for x in matrix_bigwigs(wildcards)]
    threads: 24
    conda:
        '../Envs/deeptools.yaml'
    shell:
        'computeMatrix reference-point -S {input.bigwigs} -R {input.genes} --beforeRegionStartLength 3000 --afterRegionStartLength 3000 -out {output} -p {threads} --samplesLabel {params.labels} --sortRegions keep'

rule compute_matrix_ctcf:
    input:
        bigwigs = matrix_bigwigs,
        genes = 'Macs2/CTCF_{tissue}_{rep}_summits.bed'
    output:
        'DeepTools/{tissue}_{rep}_matrix_ctcf.mat.gz'
    params:
        labels = lambda wildcards: [x.split("/")[1].split("_")[0] for x in matrix_bigwigs(wildcards)]
    threads: 24
    conda:
        '../Envs/deeptools.yaml'
    shell:
        'computeMatrix reference-point -S {input.bigwigs} -R {input.genes} --beforeRegionStartLength 500 --afterRegionStartLength 500 -out {output} -p {threads} --samplesLabel {params.labels}'

rule plot_heatmap:
    input:
        'DeepTools/{tissue}_{rep}_matrix.mat.gz'
    output:
        png = 'Figures/{tissue}_{rep}_Gene_Heatmap.png',
        clusters = 'DeepTools/{tissue}_{rep}_Gene_Clusters.bed'
    conda:
        '../Envs/deeptools.yaml'
    shell:
        'plotHeatmap -m {input} -out {output.png} --hclust 4 --colorMap plasma --outFileSortedRegions {output.clusters}'

rule plot_heatmap_ref:
    input:
        'DeepTools/{tissue}_{rep}_matrix_ref.mat.gz'
    output:
        png = 'Figures/{tissue}_{rep}_TSS_Heatmap.png',
        clusters = 'DeepTools/{tissue}_{rep}_TSS_Clusters.bed'
    threads: 12
    conda:
        '../Envs/deeptools.yaml'
    shell:
        'plotHeatmap -m {input} -out {output.png} --hclust 4 --colorMap plasma --outFileSortedRegions {output.clusters}'

rule plot_heatmap_clust:
    input:
        'DeepTools/{tissue}_{rep}_matrix_ref.mat.gz'
    output:
        png = 'Figures/{tissue}_{rep}_TSS_{clusters}_Heatmap.png',
        clusters = 'DeepTools/{tissue}_{rep}_TSS_{clusters}_Clusters.bed'
    params:
        clusters = lambda wildcards: wildcards.clusters
    threads: 12
    conda:
        '../Envs/deeptools.yaml'
    shell:
        'plotHeatmap -m {input} -out {output.png} --kmeans {params.clusters} --colorMap plasma --outFileSortedRegions {output.clusters}'
        
rule plot_heatmap_annotation:
    input:
        rules.compute_matrix_annotation.output
    output:
        png = 'Figures/{tissue}_{rep}_Annotation_Heatmap.png'
    threads: 12
    conda:
        '../Envs/deeptools.yaml'
    shell:
        'plotHeatmap -m {input} -out {output.png} --sortRegions keep --colorMap plasma --whatToShow "heatmap and colorbar"'

rule plot_heatmap_ctcf:
    input:
        'DeepTools/{tissue}_{rep}_matrix_ctcf.mat.gz'
    output:
        png = 'Figures/{tissue}_{rep}_CTCF_{clusters}_Clusters_Heatmap.png',
        clusters = 'DeepTools/{tissue}_{rep}_CTCF_{clusters}_Clusters.bed'
    params:
        clusters = lambda wildcards: wildcards.clusters
    threads: 12
    conda:
        '../Envs/deeptools.yaml'
    shell:
        'plotHeatmap -m {input} -out {output.png} --kmeans {params.clusters} --colorMap plasma --outFileSortedRegions {output.clusters}'

rule cluster_expression:
    input:
        reference = config['annotation'],
        merged = 'StringTie/Merged.gtf',
        expression = 'StringTie/RNASeq_{tissue}_{rep}_Expression.tab',
        clusters = 'DeepTools/{tissue}_{rep}_TSS_{clusters}_Clusters.bed'
    output:
        txt = 'Tables/{tissue}_{rep}_{clusters}_Clusters.txt'
    script:
        '../Scripts/Heatmap_Clusters_Expression.py'

