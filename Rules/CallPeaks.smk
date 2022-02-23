
#ruleorder: rename_atac_peaks > call_peaks

def peak_call_inputs(wildcards):
    assay, tissue, rep = wildcards.library.split('_')
    #if assay == 'DNaseSeq' or assay == 'ATAC':
    suffix = 'bam'
    #else:
    #suffix = 'tagAlign.gz'
    if rep == 'Merged':
        inputs = {'chip': ['Aligned_Reads/{}.{}'.format(x, suffix) for x in libraries(assay=assay, tissue=tissue)]}
    else:
        inputs = {'chip': ['Aligned_Reads/{}.{}'.format(wildcards.library, suffix)]}
    if assay in INPUT_ASSAYS:
        inputs['control'] = [control_library(x.split('/')[1].split('.')[0], format=suffix) for x in inputs['chip']]
    return inputs

def frip_inputs(wildcards):
    inputs = {'bam': [x.replace('.tagAlign.gz', '.bam') for x in peak_call_inputs(wildcards)['chip']]}
    inputs['bai'] = [x + '.bai' for x in inputs['bam']]
    inputs['peaks'] = 'Peak_Calls/{library}_Peaks.bed'.format(library=wildcards.library)
    return inputs

def peak_file(wildcards):
    return 'Macs2/{library}_peaks.{type}Peak'.format(library=wildcards.library, type=peak_type(wildcards.library))

def combined_peak_file(wildcards):
    return 'Peak_Calls/{assay}_{tissue}_Combined_Peaks.bed'.format(assay=wildcards.assay, tissue=wildcards.tissue)

def get_scaling_parameter(wildcards):
    import gzip
    inputs = peak_call_inputs(wildcards)
    with gzip.open(inputs['chip']) as f:
        chip_depth = 0
        for line in f:
            chip_depth += 1
    with gzip.open(inputs['control']) as f:
        input_depth = 0
        for line in f:
            input_depth += 1
    if chip_depth > input_depth:
        return '--to-large'
    else:
        return ''

##############
# Call Peaks #
##############
#ruleorder: rename_atac_peaks > call_peaks

rule call_peaks:
    input: 
        unpack(peak_call_inputs)
    output: 
        peaks = 'Peak_Calls/{library}_Peaks.bed' 
    conda:
        '../Envs/macs.yaml'
    threads: 8
    script: 
        '../Scripts/CallPeaks.py'

#rule call_atac_peaks:
    #input: 
        #chip = 'Aligned_Reads/ATAC_{sample}.bam'
    #output:
        #'Macs2/ATAC_{sample}_peaks.broadPeak' 
    #conda:
        #'../Envs/macs.yaml'
    #threads: 8
    #shell:
        #'macs2 callpeak -t {input.chip} -f BAM -n Macs2/ATAC_{wildcards.sample} -g {config[genomesize]} -q 0.05 --broad --nomodel --shift -100 --extsize 200 -B --SPMR'

#rule rename_atac_peaks:
    #input:
        #'Macs2/ATAC_{sample}_peaks.broadPeak'
    #output:
        #'Peak_Calls/ATAC_{sample}_Peaks.bed'
    #shell:
        #'cp {input} {output}'

rule fold_enrichment:
    input: 
        'Peak_Calls/{library}_Peaks.bed'
    output: 
        'Macs2/{library}_FoldEnrichment.bdg'
    conda:
        '../Envs/macs.yaml'
    threads: 8
    shell:
        'macs2 bdgcmp -t Macs2/{wildcards.library}_treat_pileup.bdg -c Macs2/{wildcards.library}_control_lambda.bdg -o {output} -m FE -p 0.000001 && '
        'sort -k1,1 -k2,2n -o {output} {output}'

rule log_likelihood:
    input: 
        'Peak_Calls/{library}_Peaks.bed'
    output: 
        'Macs2/{library}_LogLR.bdg'
    conda:
        '../Envs/macs.yaml'
    threads: 8
    shell:
        'macs2 bdgcmp -t Macs2/{wildcards.library}_treat_pileup.bdg -c Macs2/{wildcards.library}_control_lambda.bdg -o {output} -m logLR -p 0.000001 && '
        'sort -k1,1 -k2,2n -o {output} {output}'

rule overlap_peaks_with_replicate_fe:
    input:
        peakfile='Peak_Calls/{library}_Peaks.bed',
        foldenrichment=lambda wildcards: 'Macs2/{}_FoldEnrichment.bdg'.format(other_replicate(wildcards.library)),
        loglr=lambda wildcards: 'Macs2/{}_LogLR.bdg'.format(other_replicate(wildcards.library))
    output:
        'Macs2/{library}_Peak_Regions_With_Replicate_FE.txt'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        'if [ -s {input.peakfile} ]; then bedtools intersect -a {input.peakfile} -b {input.foldenrichment} -wo -sorted > {output}; else touch {output}; fi'
        
rule validate_peaks:
    input: peaks='Macs2/{library}_Peak_Regions_With_Replicate_FE.txt'
    output: outfile='Macs2/{library}_Peaks_Validated_by_Replicate.bed'
    params: peaktype=lambda wildcards: peak_type(wildcards.library)
    threads: 8
    script: '../Scripts/Validate_Peaks.py'

rule combine_peaks:
    input:
        lambda wildcards: expand('Macs2/{library}_Peaks_Validated_by_Replicate.bed', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue))
    output:
        'Peak_Calls/{assay}_{tissue}_Combined_Peaks.bed'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        'cat {input} | sort -k1,1 -k2,2n > Peak_Calls/{wildcards.assay}_{wildcards.tissue}_temp &&'
        'if [ -s Peak_Calls/{wildcards.assay}_{wildcards.tissue}_temp ]; then bedtools merge -i Peak_Calls/{wildcards.assay}_{wildcards.tissue}_temp -c 4,5 -o max > {output}; else touch {output}; fi &&'
        'rm Peak_Calls/{wildcards.assay}_{wildcards.tissue}_temp'

rule individual_frip:
    input:
        unpack(frip_inputs)
    output:
        metrics = 'Metrics/{library}_FRiP.txt',
        figure = temp('Metrics/{library}_FRiP.png')
    conda:
        '../Envs/deeptools.yaml'
    threads: 8
    script:
        '../Scripts/Calculate_FRiP.py'
        
rule combined_frip:
    input: 
        bam = lambda wildcards: expand('Aligned_Reads/{library}.bam', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue)), 
        bais = lambda wildcards: expand('Aligned_Reads/{library}.bam.bai', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue)),
        peaks = combined_peak_file
    output: metrics='Metrics/{assay}_{tissue}_FRiP.txt', figure=temp('Metrics/{assay}_{tissue}_FRiP.png')
    conda:
        '../Envs/deeptools.yaml'
    threads: 8
    shell: 'plotEnrichment --bamfiles {input.bam} --BED {input.peaks} -o {output.figure} --outRawCounts {output.metrics} -e 200 -p {threads}'
