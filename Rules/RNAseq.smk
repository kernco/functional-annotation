
#######################
# Quantify Expression #
#######################
rule stringtie:
    input: 'Aligned_Reads/{library}.bam'
    output: gtf='StringTie/{library}.gtf', tab='StringTie/{library}_Prelim_Expression.tab'
    conda: '../Envs/stringtie.yaml'
    shell: 'stringtie {input} -o {output.gtf} -p {threads} -G {config[annotation]} -l {wildcards.library} -A {output.tab}'

rule stringtie_merge:
    input: lambda wildcards: expand('StringTie/{library}.gtf', library=libraries(assay='RNASeq'))
    output: 'StringTie/Merged.gtf'
    conda: '../Envs/stringtie.yaml'
    shell: 'stringtie {input} --merge -G {config[annotation]} -o {output}'

rule stringtie_quant:
    input: gtf=rules.stringtie_merge.output, bam=rules.filter_star.output
    output: gtf='StringTie/{library}_Merged.gtf', tab='StringTie/{library}_Expression.tab'
    conda: '../Envs/stringtie.yaml'
    shell: 'stringtie {input.bam} -p {threads} -G {input.gtf} -e -o {output.gtf} -A {output.tab}'

rule stingtie_quant_annotation:
    input:
        gtf = config["annotation"],
        bam = rules.filter_star.output
    output:
        tab = 'StringTie/{library}_Annotation_Expression.tab'
    conda:
        '../Envs/stringtie.yaml'
    shell:
        'stringtie {input.bam} -p {threads} --rf -G {input.gtf} -e -A {output.tab}'

rule expression_table:
    input:
        tabs = expand('StringTie/{library}_Expression.tab', library=libraries(assay='RNASeq')),
        gtf = 'StringTie/Merged.gtf'
    output:
        tsv = 'Tables/Expression.tsv'
    script:
        '../Scripts/Make_Expression_Table.py'

rule make_tss_gtf:
    input:
        gtf = 'StringTie/Merged.gtf'
    output:
        tss = 'StringTie/Merged_TSS.bed'
    script:
        '../Scripts/Make_TSS_Gtf.py'

