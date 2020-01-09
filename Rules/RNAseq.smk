
#######################
# Quantify Expression #
#######################
#rule stringtie:
    #input: 
        #'Aligned_Reads/{library}.bam'
    #output: 
        #gtf = 'StringTie/{library}.gtf', 
        #tab = 'StringTie/{library}_Prelim_Expression.tab'
    #conda: 
        #'../Envs/stringtie.yaml'
    #shell: 
        #'stringtie {input} -o {output.gtf} -p {threads} -G {config[annotation]} -l {wildcards.library} -A {output.tab}'

#rule stringtie_merge:
    #input: 
        #lambda wildcards: expand('StringTie/{library}.gtf', library=libraries(assay='RNASeq'))
    #output: 
        #'StringTie/Merged.gtf'
    #conda: 
        #'../Envs/stringtie.yaml'
    #shell: 
        #'stringtie {input} --merge -G {config[annotation]} -o {output}'

#rule stringtie_quant:
    #input: 
        #gtf = rules.stringtie_merge.output, 
        #bam = rules.filter_star.output
    #output: 
        #gtf = 'StringTie/{library}_Merged.gtf', 
        #tab = 'StringTie/{library}_Expression.tab'
    #conda: 
        #'../Envs/stringtie.yaml'
    #shell: 
        #'stringtie {input.bam} -p {threads} -G {input.gtf} -e -o {output.gtf} -A {output.tab}'

#rule stingtie_quant_annotation:
    #input:
        #gtf = config["annotation"],
        #bam = 'Aligned_Reads/{library}.bam'
    #output:
        #tab = 'StringTie/{library}_Annotation_Expression.tab'
    #conda:
        #'../Envs/stringtie.yaml'
    #shell:
        #'stringtie {input.bam} -p {threads} --rf -G {input.gtf} -e -A {output.tab}'

#rule expression_table:
    #input:
        #tabs = expand('StringTie/{library}_Expression.tab', library=libraries(assay='RNASeq')),
        #gtf = 'StringTie/Merged.gtf'
    #output:
        #tsv = 'Tables/Novel_Expression.tsv'
    #script:
        #'../Scripts/Make_Expression_Table.py'

rule annotation_expression_table:
    input:
        tabs = lambda wildcards: expand('StringTie/{library}_Annotation_Expression.tab', library=libraries(assay='RNASeq', rep=wildcards.rep)),
        gtf = config['annotation']
    output:
        tsv = 'Tables/Annotation_Expression_{rep}.tsv'
    script:
        '../Scripts/Make_Expression_Table.py'

#rule normalize_tpms_using_mrn:
    #input:
        #tpmtables = expand('Tables/Annotation_Expression_{rep}.tsv', rep=replicates_for_assay('RNASeq'))
    #output:
        #tpmmrns = expand('Tables/Annotation_Expression_MRN_{rep}.tsv', rep=replicates_for_assay('RNASeq'))
    #script:
        #'../Scripts/MRN_Normalize_TPMs.py'

rule htseq_count:
    input:
        'Aligned_Reads/RNASeq_{sample}.bam'
    output:
        'HTSeq/RNASeq_{sample}.txt'
    conda:
        '../Envs/edger.yaml'
    threads: 8
    shell:
        'htseq-count -f bam -r pos --stranded=reverse {input} {config[annotation]} > {output}'

rule edger_normalize:
    input:
        htseqs = expand('HTSeq/{library}.txt', library=libraries(assay="RNASeq"))
    output:
        expand('Tables/Annotation_Expression_EdgeR_{rep}.tsv', rep=replicates_for_assay("RNASeq"))
    params:
        groups = [x.split('_')[1] for x in libraries(assay="RNASeq")],
        labels = ['_'.join(x.split('_')[1:]) for x in libraries(assay="RNASeq")],
        replicates = replicates_for_assay("RNASeq")
    conda:
        '../Envs/edger.yaml'
    script:
        '../Scripts/EdgeR_Normalize.R'

rule expression_pca:
    input:
        reps = expand('Tables/Annotation_Expression_EdgeR_{rep}.tsv', rep=replicates_for_assay('RNASeq'))
    output:
        png = 'Figures/Annotation_Expression_PCA.png'
    conda:
        '../Envs/r.yaml'
    script:
        '../Scripts/RNASeq_PCA.R'

rule tsi_table:
    input:
        table = 'Tables/Annotation_Expression_EdgeR_{rep}.tsv'
    output:
        table = 'Tables/Annotation_TSI_EdgeR_{rep}.tsv'
    script:
        '../Scripts/TSI_Table.py'

rule gene_expression_report:
    input:
        tpms = expand("Tables/Annotation_Expression_EdgeR_{rep}.tsv", rep=replicates_for_assay("RNASeq")),
        annotation = config['annotation']
    output:
        csv = "Tables/Gene_Expression_EdgeR_Report.csv"
    script:
        '../Scripts/Gene_Expression_Report.py'

#rule tissue_enriched_genes:
    #input:
        #reps = expand("Tables/Annotation_Expression_EdgeR_{rep}.tsv", rep=replicates_for_assay("RNASeq"))
    #output:
        #summary = "Tables/Tissue_Enriched_Genes_Summary.txt",
        #outdir = directory("Tissue_Enriched_Genes")
    #script:
        #'../Scripts/Tissue_Enriched_Genes.py'
        
rule tissue_enriched_genes:
    input:
        counts = 'Tables/Annotation_Expression_EdgeR_Average.tsv'
    output:
        outdir = directory("Tissue_Enriched_Genes")
    script:
        '../Scripts/Tissue_Enriched_Genes.py'
        
rule tpm_density_plot:
    input:
        tpms = expand("Tables/Annotation_Expression_EdgeR_{rep}.tsv", rep=replicates_for_assay("RNASeq"))
    output:
        png = "Figures/TPM_Density_Plot.png"
    conda:
        '../Envs/seaborn.yaml'
    script:
        '../Scripts/TPM_Density_Plot.py'

rule tsi_density_plot:
    input:
        tsis = expand("Tables/Annotation_TSI_EdgeR_{rep}.tsv", rep=replicates_for_assay("RNASeq"))
    output:
        png = "Figures/TSI_Density_Plot.png"
    conda:
        '../Envs/seaborn.yaml'
    script:
        '../Scripts/TSI_Density_Plot.py'

rule make_tss_gtf:
    input:
        gtf = 'StringTie/Merged.gtf'
    output:
        tss = 'StringTie/Merged_TSS.bed'
    script:
        '../Scripts/Make_TSS_Gtf.py'

rule tsi_table_mrn:
    input:
        table = 'Tables/Annotation_Expression_MRN_{rep}.tsv'
    output:
        table = 'Tables/Annotation_TSI_MRN_{rep}.tsv'
    script:
        '../Scripts/TSI_Table.py'

#rule gene_expression_report_mrn:
    #input:
        #tpms = expand("Tables/Annotation_Expression_MRN_{rep}.tsv", rep=replicates_for_assay("RNASeq")),
        #tsis = expand("Tables/Annotation_TSI_MRN_{rep}.tsv", rep=replicates_for_assay("RNASeq"))
    #output:
        #csv = "Tables/Gene_Expression_Report_MRN.csv"
    #script:
        #'../Scripts/Gene_Expression_Report.py'

#rule tpm_density_plot_mrn:
    #input:
        #tpms = expand("Tables/Annotation_Expression_MRN_{rep}.tsv", rep=replicates_for_assay("RNASeq"))
    #output:
        #png = "Figures/TPM_Density_Plot_MRN.png"
    #conda:
        #'../Envs/seaborn.yaml'
    #script:
        #'../Scripts/TPM_Density_Plot.py'

rule tsi_density_plot_mrn:
    input:
        tsis = expand("Tables/Annotation_TSI_MRN_{rep}.tsv", rep=replicates_for_assay("RNASeq"))
    output:
        png = "Figures/TSI_Density_Plot_MRN.png"
    conda:
        '../Envs/seaborn.yaml'
    script:
        '../Scripts/TSI_Density_Plot.py'
