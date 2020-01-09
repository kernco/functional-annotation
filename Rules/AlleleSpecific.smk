

#rule generate_vcf:
    #input:
        #bam = 'Aligned_Reads/{library}.bam',
        #peaks = 'Peak_Calls/{library}_Peaks.bed'
    #output:
        #'Allele_Specific/{library}.bcf'
    #conda:
        #'../Envs/samtools.yaml'
    #threads: 12
    #shell:
        #'samtools mpileup --skip-indels --positions {input.peaks} -f {config[genome]} -g -o {output} {input.bam}'

#rule call_snps:
    #input:
        #'Allele_Specific/{library}.bcf'
    #output:
        #'Allele_Specific/{library}_SNPs'
    #conda:
        #'../Envs/bcftools.yaml'
    #threads: 12
    #shell:
        #'bcftools call -mv --threads {threads} -A -o {output} {input}'

def input_libs(wildcards):
    for assay in config['inputs']:
        for library in libraries(assay=assay, rep=wildcards.rep, skip_input=False):
            yield 'Aligned_Reads/{library}.bam'.format(library=library)

rule make_genome_dict:
    input:
        genome = config["genome"]
    output:
        config["genome"].replace('.fa', '.dict')
    conda:
        '../Envs/picard.yaml'
    shell:
        'picard-tools CreateSequenceDictionary R= {input.genome} O= {output}'

rule combine_input:
    input:
        bams = input_libs
    output:
        temp('Aligned_Reads/AllInput_{rep}_temp.bam')
    conda:
        '../Envs/samtools.yaml'
    shell:
        'samtools merge {output} {input.bams}'

rule combine_read_groups:
    input:
        rules.combine_input.output
    output:
        'Aligned_Reads/AllInput_{rep}.bam'
    conda:
        '../Envs/picard.yaml'
    shell:
        'picard-tools AddOrReplaceReadGroups I={input} O={output} RGLB=AllInput_{wildcards.rep} RGPL=illumina RGPU=run1 RGSM=AllInput_{wildcards.rep}'

rule gatk_haplotype:
    input:
        ref = config["genome"],
        dct = config["genome"].replace('.fa', '.dict'),
        bam = 'Aligned_Reads/{name}.bam',
        bai = 'Aligned_Reads/{name}.bam.bai'
    output:
        vcf = 'Allele_Specific/{name}.vcf'
    conda:
        '../Envs/gatk.yaml'
    threads: 8
    shell:
        'gatk --java-options "-Xmx64G" HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf}'

rule gatk_filter:
    input:
        vcf = 'Allele_Specific/{name}.vcf',
        ref = config["genome"]
    output:
        vcf = 'Allele_Specific/{name}.filtered.vcf'
    conda:
        '../Envs/gatk.yaml'
    threads: 8
    shell:
        'gatk --java-options "-Xmx8G" VariantFiltration -R {input.ref} -V {input.vcf} -O {output.vcf} --filter-expression "QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -filter-name "filter1"'

rule overlap_snps_with_peaks:
    input:
        vcf = 'Allele_Specific/{assay}_{tissue}_{rep}.filtered.vcf',
        peaks = 'Peak_Calls/{assay}_{tissue}_Combined_Peaks.bed'
    output:
        'Allele_Specific/{assay}_{tissue}_{rep}.Peak_SNPs.vcf'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        'bedtools intersect -a {input.vcf} -b {input.peaks} -u > {output}'

rule create_baal_sample_sheet:
    output:
        tsv = 'Allele_Specific/Sample_Sheet_{assay}.tsv'
    run:
        with open(output.tsv, 'w') as f:
            print('group_name\ttarget\treplicate_number\tbam_name\tbed_name', file=f)
            for library in libraries(assay=wildcards.assay):
                assay, tissue, rep = library.split('_')
                f.write('GRP\t{}_{}\t{}\t{}\t{}\n'.format(tissue, assay, rep, 'Aligned_Reads/{}.bam'.format(library), 'Peak_Calls/{}_{}_Combined_Peaks.bed'.format(assay, tissue)))

rule filter_snps:
    input:
        bed = 'Allele_Specific/{file}.bed',
        peaks = lambda wildcards: expand('Peak_Calls/{assay}_{tissue}_Combined_Peaks.bed', assay=wildcards.assay, tissue=config['tissues'])
    output:
        'Allele_Specific/{file}_{assay}.bed'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        'bedtools intersect -a {input.bed} -b {input.peaks} > {output}'
    
rule bed_to_hets:
    input:
        vcf = 'Allele_Specific/{file}_{assay}.bed'
    output:
        hets = 'Allele_Specific/{file}_{assay}.hets'
    script:
        '../Scripts/bed2hets.py'

rule allele_specific_binding:
    input:
        samplesheet = 'Allele_Specific/Sample_Sheet_{assay}.tsv',
        hets = 'Allele_Specific/SNPs_{assay}.hets',
        gdna = libraries(assay='Input')
    output:
        table = 'Allele_Specific/ASB_{assay}.txt',
        qc = 'Allele_Specific/Filtering_QC_{assay}.txt',
        summary = 'Allele_Specific/Summary_{assay}.txt'
    conda:
        '../Envs/r.yaml'
    threads: 24
    script:
        '../Scripts/AlleleSpecificBinding.R'
