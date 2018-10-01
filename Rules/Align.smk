
ruleorder: star_align > bwa_mem_paired
ruleorder: filter_star > add_read_group #remove_duplicates
ruleorder: trim_atac_paired > trim_paired_reads

def alignment_report_inputs(wildcards):
    if wildcards.assay.startswith("RNASeq"):
        return expand(['Aligned_Reads/{library}.bam', 'Trimmed_Reads/{library}_R2.fq.gz_trimming_report.txt', 'Trimmed_Reads/{library}_R2_val_2.fq.gz', 'STAR_Output/{library}_Log.final.out'], library=libraries(assay=wildcards.assay))
    else:
        return expand(['Aligned_Reads/{library}.bam', 'Bwa_Output/{library}.aligned.bam', 'Bwa_Output/{library}.filtered.bam', 'Bwa_Output/{library}.duplicate-marked.bam'], library=libraries(assay=wildcards.assay, skip_input=False))

#######
# BWA #
#######
rule trim_reads:
    input: 
        'Raw_Reads/{library}.fq.gz'
    output: 
        'Trimmed_Reads/{library}_trimmed.fq.gz', 
        'Trimmed_Reads/{library}.fq.gz_trimming_report.txt'
    conda:
        '../Envs/trimgalore.yaml'
    shell: 
        'trim_galore -q 20 {input} -o Trimmed_Reads'

rule trim_paired_reads:
    input: 
        'Raw_Reads/{library}_R1.fq.gz', 
        'Raw_Reads/{library}_R2.fq.gz'
    output: 
        'Trimmed_Reads/{library}_R1_val_1.fq.gz', 
        'Trimmed_Reads/{library}_R2_val_2.fq.gz',
        'Trimmed_Reads/{library}_R1.fq.gz_trimming_report.txt',
        'Trimmed_Reads/{library}_R2.fq.gz_trimming_report.txt'
    conda:
        '../Envs/trimgalore.yaml'
    shell: 
        'trim_galore -q 20 {input} -o Trimmed_Reads --paired'

rule trim_atac_paired:
    input:
        'Raw_Reads/ATAC_{sample}_R1.fq.gz',
        'Raw_Reads/ATAC_{sample}_R2.fq.gz'
    output:
        'Trimmed_Reads/ATAC_{sample}_R1_val_1.fq.gz',
        'Trimmed_Reads/ATAC_{sample}_R2_val_2.fq.gz'
    conda:
        '../Envs/trimgalore.yaml'
    shell:
        'trim_galore -q 20 -a CTGTCTCTTATA -stringency 1 --paired --retain_unpaired --length 10 -o Trimmed_Reads {input}'

rule bwa_index:
    input: 
        config["genome"]
    output: 
        config["genome"] + ".amb",
        config["genome"] + ".ann",
        config["genome"] + ".bwt",
        config["genome"] + ".pac",
        config["genome"] + ".sa"
    params:
        prefix = config["genome"],
        algorithm = "bwtsw"
    conda:
        '../Envs/bwa.yaml'
    #wrapper:
        #"0.23.1/bio/bwa/index"
    shell:
        'bwa index -a {params.algorithm} {input}'

rule bwa_mem:
    input: 
        rules.bwa_index.output,
        reads = ['Trimmed_Reads/{library}_trimmed.fq.gz']
    output: 
        'Bwa_Output/{library}.aligned.bam'
    params:
        index = config["genome"],
        extra = r"-M",
        sort = "samtools",
        sort_order = "coordinate"
    threads: 24
    conda:
        '../Envs/bwa.yaml'
    #wrapper:
        #"0.23.1/bio/bwa/mem"
    shell: 
        'bwa mem -M -t {threads} {config[genome]} {input.reads} | samtools view -bS - > {output}'

rule bwa_mem_paired:
    input: 
        rules.bwa_index.output, 
        reads = ['Trimmed_Reads/{library}_R1_val_1.fq.gz', 'Trimmed_Reads/{library}_R2_val_2.fq.gz']
    output: 
        'Bwa_Output/{library}.aligned.bam'
    params:
        index = config["genome"],
        extra = r"-M",
        sort = "samtools",
        sort_order = "coordinate"
    threads: 24
    conda:
        '../Envs/bwa.yaml'
    #wrapper:
        #"0.23.1/bio/bwa/mem"
    shell: 
        'bwa mem -M -t {threads} {config[genome]} {input.reads} | samtools view -bS - > {output}'


rule filter_alignments:
    input: 
        'Bwa_Output/{library}.aligned.bam'
    output: 
        'Bwa_Output/{library}.filtered.bam'
    conda:
        '../Envs/samtools.yaml'
    shell: 
        'samtools view -h -F 1804 -q {config[mapq]} {input} | grep -v XA:Z | grep -v SA:Z | samtools view -S -b - > {output}'

rule sort_bam:
    input: 
        rules.filter_alignments.output
    output: 
        'Bwa_Output/{library}.sorted.bam'
    threads: 12
    conda:
        '../Envs/samtools.yaml'
    shell: 
        'samtools sort -@ {threads} -o Bwa_Output/{wildcards.library}.sorted.bam {input}'

rule mark_duplicates:
    input: 
        rules.sort_bam.output
    output: 
        bam = 'Bwa_Output/{library}.duplicate-marked.bam', 
        metrics = 'Metrics/{library}.dup_metrics'
    log:
        'Logs/Picard/Dedup/{library}.log'
    #params:
        #"REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT"
    #wrapper:
        #"0.23.1/bio/picard/markduplicates"
    threads: 24
    conda:
        '../Envs/picard.yaml'
    shell: 
        'picard-tools MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT'

rule remove_duplicates:
    input: 
        rules.mark_duplicates.output.bam
    output: 
        #'Bwa_Output/{library}.no-duplicates.bam',
        #bam = 'Aligned_Reads/{library}.bam'
        bam = temp('Bwa_Output/{library}.no-duplicates.bam')
    conda:
        '../Envs/samtools.yaml'
    shell: 
        'samtools view -F 1804 -q {config[mapq]} -b {input} > {output.bam}'

rule add_read_group:
    input:
        rules.remove_duplicates.output.bam
    output:
        'Aligned_Reads/{library}.bam'
    conda:
        '../Envs/picard.yaml'
    threads: 12
    shell:
        'picard-tools AddOrReplaceReadGroups I={input} O={output} RGLB={wildcards.library} RGPL=illumina RGPU=l1 RGSM={wildcards.library}'

#rule remove_mitochondrial:
    #input:
        #bam = 'Bwa_Output/{library}.duplicate-marked.bam',
        #idx = 'Bwa_Output/{library}.duplicate-marked.bam.bai'
    #output:
        #'Aligned_Reads/{library}.bam'
    #conda:
        #'../Envs/samtools.yaml'
    #shell:
        #'samtools idxstats {input.bam} | cut -f 1 | grep -v MT | xargs samtools view -b {input.bam} > {output}'

rule bam_to_tagalign:
    input: 
        'Aligned_Reads/{library}.bam'
    output: 
        'Aligned_Reads/{library}.tagAlign.gz'
    conda:
        '../Envs/bedtools.yaml'
    threads: 8
    shell: 
        """bedtools bamtobed -i {input} | awk 'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}' | shuf -n50000000 | gzip -nc > {output}"""
    
rule index_bam:
    input: 
        'Aligned_Reads/{name}.bam'
    output: 
        'Aligned_Reads/{name}.bam.bai'
    #wrapper:
        #"0.23.1/bio/samtools/index"
    conda:
        '../Envs/samtools.yaml'
    shell: 
        'samtools index {input}'

########
# STAR #
########
rule star_index:
    input: 
        genome = config['genome'], 
        gtf = config['annotation']
    output: 
        'Genome_Index/chrName.txt'
    conda:
        '../Envs/star.yaml'
    threads: 24
    shell: 
        'mkdir -p Genome_Index && STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir Genome_Index --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf}'

rule star_align:
    input: 
        r1 = 'Trimmed_Reads/{library}_R1_val_1.fq.gz', 
        r2 = 'Trimmed_Reads/{library}_R2_val_2.fq.gz', 
        index = 'Genome_Index/chrName.txt'
    output: 
        bam='STAR_Output/{library}_Aligned.sortedByCoord.out.bam',
        log='STAR_Output/{library}_Log.final.out'
    wildcard_constraints: 
        library="RNASeq_.+"
    conda:
        '../Envs/star.yaml'
    threads: 24
    shell: 
        'STAR --runThreadN {threads} --genomeDir Genome_Index --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat --outFileNamePrefix STAR_Output/{wildcards.library}_ --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --outSAMtype BAM SortedByCoordinate'

rule filter_star:
    input: 
        rules.star_align.output.bam
    output: 
        'Aligned_Reads/{library}.bam'
    wildcard_constraints: library="RNASeq_.+"
    conda:
        '../Envs/samtools.yaml'
    shell: 
        'samtools view -b -q {config[mapq]} {input} > {output}'

rule make_alignment_report:
    input: 
        alignment_report_inputs
    output: 
        txt = 'Tables/{assay}_Alignment_Summary.txt', 
        csv = 'Tables/{assay}_Alignment_Summary.csv'
    params: 
        libraries = lambda wildcards: expand('{library}', library=libraries(assay=wildcards.assay, skip_input=False))
    script: 
        '../Scripts/Make_Alignment_Summary.py'
