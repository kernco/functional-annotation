rule trim_methylation:
    input:
        'Raw_Reads/RRBS_{sample}.fq.gz'
    output:
        'Trimmed_Reads/RRBS_{sample}_trimmed.fq.gz',
        'Trimmed_Reads/RRBS_{sample}.fq.gz_trimming_report.txt'
    shell:
        'module load trim-galore python cutadapt && trim_galore --rrbs {input} -o Trimmed_Reads'

rule bismark_prepare_genome:
    input:
        config["genome"]
    output:
        'Genome/Bisulfite_Genome'
    shell:
        'bismark_genome_preparation --genomic_composition Genome'

rule bismark:
    input:
        'Trimmed_Reads/RRBS_{sample}_trimmed.fq.gz',
        'Genome/Bisulfite_Genome'
    output:
        'Bismark/RRBS_{sample}_trimmed.fq.gz_bismark.bam',
        'Bismark/RRBS_{sample}_trimmed.fq.gz_bismark_SE_report.txt'
    threads: 4
    shell:
        'bismark --genome Genome --multicore {threads} -o Bismark {input}'
