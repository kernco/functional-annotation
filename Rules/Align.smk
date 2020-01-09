
##
# Trim reads 
##
rule trim_reads:
    input: 
        'Raw_Reads/{library}.fq.gz'
    output: 
        reads = temp(config['tempdir'] + '/{library}_trimmed.fq.gz'),
        report = temp(config['tempdir'] + '/{library}.fq.gz_trimming_report.txt')
    #group: 'bwa_align'
    conda:
        '../Envs/trimgalore.yaml'
    shell: 
        'trim_galore -q 20 {input} -o {config[tempdir]}'

rule trim_paired_reads:
    input: 
        'Raw_Reads/{library}_R1.fq.gz', 
        'Raw_Reads/{library}_R2.fq.gz'
    output: 
        reads = temp(config['tempdir'] + '/{library}_R1_val_1.fq.gz'), 
        reads2 = temp(config['tempdir'] + '/{library}_R2_val_2.fq.gz'),
        report = temp(config['tempdir'] + '/{library}_R1.fq.gz_trimming_report.txt'),
        report2 = temp(config['tempdir'] + '/{library}_R2.fq.gz_trimming_report.txt')
    params:
        trim_params = lambda wildcards: '-a CTGTCTCTTATA -stringency 1 --retain_unpaired --length 10' if wildcards.library.startswith('ATAC') else ''
    #group: 'bwa_align'
    conda:
        '../Envs/trimgalore.yaml'
    shell: 
        'trim_galore -q 20 --paired {params.trim_params} {input} -o {config[tempdir]}'

##
# Index genome
##
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
    threads: 24
    shell:
        'bwa index -a {params.algorithm} {input}'

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

##
# Align reads
#
# Single-end workflow: bwa_mem
# Paired-end workflow: bwa_mem_paired
# RNA-specific workflow: star_align
##
#
ruleorder: star_align >  bwa_mem_paired
ruleorder: star_align_single > bwa_mem
ruleorder: star_align_single > star_align
ruleorder: bwa_mem > bwa_mem_paired

localrules: index_bam
rule index_bam:
    input: 
        '{name}.bam'
    output: 
        '{name}.bam.bai'
    #conda:
        #'../Envs/samtools.yaml'
    shell: 
        'samtools index {input}'

rule bwa_mem:
    input: 
        index = rules.bwa_index.output,
        reads = rules.trim_reads.output.reads
    output: 
        bam = temp(config['tempdir'] + '/{library}.aligned.bam')
    #group: 'bwa_align'
    threads: 24
    conda:
        '../Envs/bwa.yaml'
    shell: 
        "bwa mem -M -t {threads} -R '@RG\\tID:{wildcards.library}\\tSM:{wildcards.library}' {config[genome]} {input.reads} | samtools view -bS - > {output}"

#Temporary to make RNA-seq work
rule bwa_mem_paired:
    input: 
        rules.bwa_index.output, 
        reads = [config['tempdir'] + '/{library}_R1_val_1.fq.gz', config['tempdir'] + '/{library}_R2_val_2.fq.gz']
    output: 
        bam = temp(config['tempdir'] + '/{library}.aligned.bam')
    #group: 'bwa_align'
    threads: 24
    conda:
        '../Envs/bwa.yaml'
    shell: 
        "bwa mem -M -t {threads} -R '@RG\\tID:{wildcards.library}\\tSM:{wildcards.library}' {config[genome]} {input.reads} | samtools view -bS - > {output}"

rule star_align:
    input: 
        r1 = config['tempdir'] + '/{library}_R1_val_1.fq.gz', 
        r2 = config['tempdir'] + '/{library}_R2_val_2.fq.gz', 
        index = 'Genome_Index/chrName.txt'
    output: 
        bam = temp(config['tempdir'] + '/{library}_Aligned.sortedByCoord.out.bam'),
        log = config['tempdir'] + '/{library}_Log.final.out'
    wildcard_constraints: 
        library="RNASeq_.+"
    conda:
        '../Envs/star.yaml'
    #group: "bwa_align"
    threads: 24
    shell: 
        'STAR --runThreadN {threads} --genomeDir Genome_Index --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat --outFileNamePrefix {config[tempdir]}/{wildcards.library}_ --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --outSAMtype BAM SortedByCoordinate'

rule star_align_single:
    input:
        r1 = rules.trim_reads.output.reads,
        index = 'Genome_Index/chrName.txt'
    output:
        bam = temp(config['tempdir'] + '/{library}_Aligned.sortedByCoord.out.bam'),
        log = config['tempdir'] + '/{library}_Log.final.out'
    wildcard_constraints:
        library="RNASeq_.+"
    conda:
        '../Envs/star.yaml'
    #group: "bwa_align"
    threads: 24
    shell: 
        'STAR --runThreadN {threads} --genomeDir Genome_Index --readFilesIn {input.r1} --readFilesCommand zcat --outFileNamePrefix {config[tempdir]}/{wildcards.library}_ --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --outSAMtype BAM SortedByCoordinate'
        
##
# Filter alignments
#
# Default workflow: filter_alignments -> sort_bam
# ATAC-specific workflow: filter_atac_alignments -> sort_bam
# RNA-specific workflow: filter_star
##
ruleorder: filter_star > remove_duplicates
ruleorder: filter_atac_alignments > filter_alignments
ruleorder: filter_atac_alignments > remove_duplicates

rule filter_alignments:
    input: 
        config['tempdir'] + '/{library}.aligned.bam'
    output: 
        bam = temp(config['tempdir'] + '/{library}.filtered.bam')
    #group: 'bwa_align'
    conda:
        '../Envs/samtools.yaml'
    shell: 
        'samtools view -h -F 1804 -q {config[mapq]} {input} | grep -v XA:Z | grep -v SA:Z | samtools view -S -b - > {output}'

rule filter_atac_alignments:
    input:
        bam = config['tempdir'] + '/{library}.aligned.bam',
    output:
        intermediate = temp(config['tempdir'] + '/{library}.intermediate.bam'),
        bam = temp(config['tempdir'] + '/{library}.filtered.bam')
    wildcard_constraints:
        library = 'ATAC_.+'
    #group: 'bwa_align'
    threads: 12
    conda:
        '../Envs/samtools.yaml'
    shell:
        'samtools sort -@ {threads} -o {input.bam}.sorted {input.bam} && samtools index {input.bam}.sorted &&'
        'samtools idxstats {input.bam}.sorted | cut -f 1 | grep -v chrM | xargs samtools view -b {input.bam}.sorted > {output.intermediate} && '
        'samtools view -h -F 1804 -q {config[mapq]} {output.intermediate} | grep -v XA:Z | grep -v SA:Z | samtools view -S -b - > {output.bam}'

rule filter_star:
    input: 
        rules.star_align.output.bam
    output: 
        temp(config['tempdir'] + '/{library}.final.bam')
    wildcard_constraints: 
        library = "RNASeq_.+"
    #group: "bwa_align"
    conda:
        '../Envs/samtools.yaml'
    shell: 
        'samtools view -b -q {config[mapq]} {input} > {output}'

rule namesort_bam:
    input:
        config['tempdir'] + '/{library}.filtered.bam'
    output:
        temp(config['tempdir'] + '/{library}.namesorted.bam')
    #group: 'bwa_align'
    threads: 24
    conda:
        '../Envs/samtools.yaml'
    shell:
        'samtools sort -@ {threads} -n -o {output} {input}'

rule fixmate_bam:
    input:
        config['tempdir'] + '/{library}.namesorted.bam'
    output:
        temp(config['tempdir'] + '/{library}.fixmate.bam')
    #group: 'bwa_align'
    threads: 24
    conda:
        '../Envs/samtools.yaml'
    shell:
        'samtools fixmate -m {input} {output}'

rule sort_bam:
    input: 
        config['tempdir'] + '/{library}.fixmate.bam'
    output: 
        temp(config['tempdir'] + '/{library}.sorted.bam')
    #group: 'bwa_align'
    threads: 24
    conda:
        '../Envs/samtools.yaml'
    shell:
        'samtools sort -@ {threads} -o {output} {input}'

##
# Mark and remove duplicates
#
# Default workflow: mark_duplicates -> remove_duplicates
# RNA-specific workflow: Nothing
##
rule mark_duplicates:
    input:
        rules.sort_bam.output
    output:
        bam = temp(config['tempdir'] + '/{library}.duplicate-marked.bam')
    #group: 'bwa_align'
    threads: 24
    conda:
        '../Envs/samtools.yaml'
    shell:
        'samtools markdup {input} {output.bam}'

#rule mark_duplicates:
    #input: 
        #rules.sort_bam.output
    #output: 
        #bam = temp(config['tempdir'] + '/{library}.duplicate-marked.bam'),
        #metrics = temp(config['tempdir'] + '/{library}.dup_metrics')
    #group: 'bwa_align'
    #log:
        #'Logs/Picard/Dedup/{library}.log'
    #threads: 24
    #conda:
        #'../Envs/picard.yaml'
    #shell: 
        #'picard-tools MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR={config[tempdir]}'

rule remove_duplicates:
    input: 
        rules.mark_duplicates.output.bam
    output: 
        bam = temp(config['tempdir'] + '/{library}.final.bam')
    #group: 'bwa_align'
    conda:
        '../Envs/samtools.yaml'
    shell: 
        'samtools view -F 1804 -q {config[mapq]} -b {input} > {output.bam}'

################################
# Create final alignment files #
################################
#rule add_read_group:
    #input:
        #rules.remove_duplicates.output.bam
    #output:
        #temp(config['tempdir'] + '/{library}.final.bam')
    #group: 'bwa_align'
    #conda:
        #'../Envs/picard.yaml'
    #threads: 12
    #shell:
        #'picard-tools AddOrReplaceReadGroups I={input} O={output} RGLB={wildcards.library} RGPL=illumina RGPU=l1 RGSM={wildcards.library} RGID={wildcards.library}'

rule bam_to_tagalign:
    input: 
        config['tempdir'] + '/{library}.final.bam'
    output: 
        temp(config['tempdir'] + '/{library}.tagAlign.gz')
    #group: 'bwa_align'
    conda:
        '../Envs/bedtools.yaml'
    threads: 8
    shell: 
        """bedtools bamtobed -i {input} | awk 'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}' | gzip -nc > {output}"""

##########################
# Calculate some metrics #
##########################
rule spp_stats:
    input: 
        rules.bam_to_tagalign.output
    output: 
        stats = 'Metrics/{library}.spp_stats.txt', 
        figure = 'Metrics/{library}_Cross_Correlation.pdf'
    #group: 'bwa_align'
    threads: 24
    conda:
        '../Envs/r.yaml'
    shell: 
        'Rscript /home/ckern/phantompeakqualtools/run_spp.R -c={input} -rf -out={output.stats} -p={threads} -s=0:2:400 -savp={output.figure} -tmpdir={config[tempdir]}'

############################
# Get alignment statistics #
############################

def alignment_stats_inputs(wildcards):
    inputs = {}
    if wildcards.library in PAIRED_LIBS:
        inputs['trimmed_fq'] = config['tempdir'] + '/{library}_R1_val_1.fq.gz'.format(library=wildcards.library)
        inputs['trim_report'] = config['tempdir'] + '/{library}_R1.fq.gz_trimming_report.txt'.format(library=wildcards.library)
    else:
        inputs['trimmed_fq'] = config['tempdir'] + '/{library}_trimmed.fq.gz'.format(library=wildcards.library)
        inputs['trim_report'] = config['tempdir'] + '/{library}.fq.gz_trimming_report.txt'.format(library=wildcards.library)
    if not wildcards.library.startswith('RNASeq'):
        inputs['aligned_bam'] = config['tempdir'] + '/{library}.aligned.bam'.format(library=wildcards.library)
        inputs['filtered_bam'] = config['tempdir'] + '/{library}.filtered.bam'.format(library=wildcards.library)
        inputs['deduped_bam'] = config['tempdir'] + '/{library}.final.bam'.format(library=wildcards.library)
    else:
        inputs['star_log'] = config['tempdir'] + '/{library}_Log.final.out'.format(library=wildcards.library)
        inputs['final_bam'] = config['tempdir'] + '/{library}.final.bam'.format(library=wildcards.library)
    if (not (wildcards.library.startswith('ATAC')) and (not wildcards.library.startswith('RNASeq'))):
        inputs['spp_stats'] = 'Metrics/{library}.spp_stats.txt'.format(library=wildcards.library)
    return inputs

def get_assay_type(wildcards):
    #TODO: Make this customizable in config file
    if wildcards.library.startswith('ATAC'):
        return 'ATAC'
    elif wildcards.library.startswith('DNase'):
        return 'Paired'
    elif wildcards.library.startswith('RNASeq'):
        return 'RNA-seq'
    else:
        return 'Single'

rule get_alignment_stats:
    input:
        unpack(alignment_stats_inputs)
    output:
        json = 'Metrics/{library}_Alignment_Stats.json'
    params:
        assay_type = get_assay_type
    #group: 'bwa_align'
    script:
        '../Scripts/Get_Alignment_Stats.py'

def align_library_inputs(wildcards):
    inputs = {}
    if not wildcards.library.startswith('ATAC'):
        inputs['stats'] = 'Metrics/{}_Alignment_Stats.json'.format(wildcards.library)
    inputs['final_bam'] = config['tempdir'] + '/{}.final.bam'.format(wildcards.library)
    if wildcards.library not in PAIRED_LIBS and not wildcards.library.startswith('RNASeq'):
        inputs['tagalign'] = config['tempdir'] + '/{}.tagAlign.gz'.format(wildcards.library)
    return inputs
    
rule align_library:
    input:
        unpack(align_library_inputs)
    output:
        final_bam = 'Aligned_Reads/{library}.bam',
        tagalign = 'Aligned_Reads/{library}.tagAlign.gz'
    #group: 'bwa_align'
    run:
        shell('cp {input.final_bam} {output.final_bam}')
        if hasattr(input, 'tagalign'):
            shell('cp {input.tagalign} {output.tagalign}')
        else:
            shell('touch {output.tagalign}')

###############################
# Summarize alignment reports #
###############################
rule make_alignment_report:
    input:
        stats = lambda wildcards: expand('Metrics/{library}_Alignment_Stats.json', library=libraries(wildcards.assay, skip_input=False)),
        bams = lambda wildcards: expand('Aligned_Reads/{library}.bam', library=libraries(wildcards.assay, skip_input=False))
    output:
        txt = 'Tables/{assay}_Alignment_Summary.txt',
        csv = 'Tables/{assay}_Alignment_Summary.csv'
    params:
        libraries = lambda wildcards: expand('{library}', library=libraries(assay=wildcards.assay, skip_input=False))
    script:
        '../Scripts/Make_Alignment_Summary.py'
