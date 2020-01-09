rule trim_methylation:
    input:
        'Raw_Reads/RRBS_{sample}.fq.gz'
    output:
        temp('Temp/RRBS_{sample}_trimmed.fq.gz'),
        temp('Temp/RRBS_{sample}.fq.gz_trimming_report.txt')
    conda:
        '../Envs/trimgalore.yaml'
    shell:
        'trim_galore --rrbs {input} -o {config[tempdir]}'

rule bismark_prepare_genome:
    input:
        config["genome"]
    output:
        directory(os.path.dirname(config["genome"]) + '/Bisulfite_Genome')
    params:
        genomedir = os.path.dirname(config["genome"])
    conda:
        '../Envs/bismark.yaml'
    threads: 12
    shell:
        'bismark_genome_preparation --genomic_composition {params.genomedir}'

rule bismark:
    input:
        genome = os.path.dirname(config["genome"]) + '/Bisulfite_Genome',
        reads = config['tempdir'] + '/RRBS_{sample}_trimmed.fq.gz'
    output:
        'Bismark/RRBS_{sample}_trimmed_bismark_bt2.bam',
        'Bismark/RRBS_{sample}_trimmed_bismark_bt2_SE_report.txt'
    params:
        genomdir = os.path.dirname(config["genome"])
    conda:
        '../Envs/bismark.yaml'
    threads: 24
    shell:
        'bismark -p 8 -o Bismark {params.genomdir} {input.reads}'

rule bs_seeker_build:
    input:
        genome = config["genome"],
    output:
        directory('BSSeeker_Genome')
    conda:
        '../Envs/bsseeker.yaml'
    threads: 24
    shell:
        'bs_seeker2-build.py --rrbs -f {input.genome} --aligner=bowtie2 -c C-CGG -d BSSeeker_Genome'

ruleorder: bs_seeker_align > align_library
rule bs_seeker_align:
    input:
        fastq = 'Temp/RRBS_{sample}_trimmed.fq.gz',
        genome = config["genome"]
    output:
        bam = 'Aligned_Reads/RRBS_{sample}.bam',
        log = 'Metrics/RRBS_{sample}_Align.log'
    conda:
        '../Envs/bsseeker.yaml'
    threads: 24
    shell:
        'bs_seeker2-align.py --rrbs -i {input.fastq} -c C-CGG,T-CGA -g {input.genome} -d BSSeeker_Genome --aligner=bowtie2 -o {output.bam} --bt2-p 24 --bt2--mm > {output.log}'
        
rule bs_seeker_call:
    input:
        bam = 'Aligned_Reads/RRBS_{sample}.bam'
    output:
        wig = 'Methylation/{sample}.wig',
        cgmap = 'Methylation/{sample}.CGmap.gz',
        atcgmap = 'Methylation/{sample}.ATCGmap.gz'
    conda:
        '../Envs/bsseeker.yaml'
    threads: 8
    shell:
        'bs_seeker2-call_methylation.py -i {input.bam} -d BSSeeker_Genome/susScr11.fa_rrbs_CCGG-TCGA_20_500_bowtie2 -o Methylation/{wildcards.sample} -x --rm-CCGG'

rule merge_methylation_calls:
    input:
        reps = lambda wildcards: expand('Methylation/{tissue}_{rep}.CGmap.gz', tissue=wildcards.tissue, rep=replicates_for_assay('RRBS'))
    output:
        'Methylation/{tissue}_Merged.CGmap.gz'
    conda:
        '../Envs/macs.yaml' #This is so it runs in a python2 environment
    shell:
        'cgmaptools merge2 cgmap -1 {input.reps[0]} -2 {input.reps[1]} -o {output}'

rule call_to_bed:
    input:
        cgmap = 'Methylation/{sample}.CGmap.gz'
    output:
        bed = 'Methylation/{sample}.bed'
    run:
        import gzip
        outfile = open(output.bed, 'w')
        with gzip.open(input.cgmap) as infile:
            for line in infile:
                cols = line.decode('utf-8').split()
                if cols[3] == 'CG' and int(cols[7]) >= 10:
                    print("\t".join([cols[0], cols[2], str(int(cols[2]) + 1), cols[5], cols[6], cols[7]]), file=outfile)

