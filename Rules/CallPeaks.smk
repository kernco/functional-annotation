
ruleorder: rename_atac_peaks > call_peaks

def peak_call_inputs(wildcards):
    inputs = {'chip': 'Aligned_Reads/{}.tagAlign.gz'.format(wildcards.library)}
    if wildcards.library.split('_')[0] not in config['no_input']:
        inputs['control'] = control_library(wildcards.library, format="tagAlign")
    return inputs

def idr_peak_call_inputs(wildcards):
    inputs = {'chip': 'IDR/{}.tagAlign.gz'.format(wildcards.idrprefix)}
    assay = wildcards.idrprefix.split('_')[0]
    if assay not in config['no_input']:
        if 'Pool' in wildcards.idrprefix:
            assay, tissue = wildcards.idrprefix.split('.')[0].split('_')
            inputs['control'] = 'IDR/{assay}_{tissue}.Pooled.tagAlign.gz'.format(assay=config[assay + "_input"], tissue=tissue)
        else:
            assay, tissue, rep = wildcards.idrprefix.split('_')
            inputs['control'] = 'IDR/{}.tagAlign.gz'.format('_'.join([config[assay + "_input"], tissue, rep]))
    return inputs

def peak_file(wildcards):
    return 'Macs2/{library}_peaks.{type}Peak'.format(library=wildcards.library, type=peak_type(wildcards.library))

def combined_peak_file(wildcards):
    return 'Peak_Calls/{assay}_{tissue}_Combined_Peaks.bed'.format(assay=wildcards.assay, tissue=wildcards.tissue)

def idr_peak_file(wildcards):
    return 'Peak_Calls/{assay}_{tissue}_IDR.{type}Peak'.format(assay=wildcards.assay, tissue=wildcards.tissue, type=peak_type(wildcards.assay))

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

#######
# IDR #
#######

rule create_pseudoreplicates:
    input: 'Aligned_Reads/{library}.tagAlign.gz'
    output: pr1='IDR/{library}.pr1.tagAlign.gz', pr2='IDR/{library}.pr2.tagAlign.gz'
    shell:
        """
        nlines=$( zcat {input} | wc -l )
        nlines=$(( (nlines + 1) / 2 ))
        zcat {input} | shuf --random-source={input} | split -d -l ${{nlines}} - IDR/{wildcards.library}
        gzip -nc "IDR/{wildcards.library}00" > {output.pr1}
        rm "IDR/{wildcards.library}00"
        gzip -nc "IDR/{wildcards.library}01" > {output.pr2}
        rm "IDR/{wildcards.library}01"
        """

rule create_pooled:
    input: 
        lambda wildcards: expand('Aligned_Reads/{library}.tagAlign.gz', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue, skip_input=False))
    output: 
        pool='IDR/{assay}_{tissue}.Pooled.tagAlign.gz'
    shell:
        'zcat {input} | gzip -nc > {output}'

rule create_pooled_pr1:
    input:
        lambda wildcards: expand('IDR/{library}.pr1.tagAlign.gz', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue, skip_input=False))
    output:
        'IDR/{assay}_{tissue}.Poolr1.tagAlign.gz'
    shell:
        'zcat {input} | gzip -nc > {output}'

rule create_pooled_pr2:
    input:
        lambda wildcards: expand('IDR/{library}.pr2.tagAlign.gz', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue, skip_input=False))
    output:
        'IDR/{assay}_{tissue}.Poolr2.tagAlign.gz'
    shell:
        'zcat {input} | gzip -nc > {output}'
        
rule call_narrow_idr_peaks:
    input:
        unpack(idr_peak_call_inputs)
    output:
        'IDR/{idrprefix}_peaks.narrowPeak'
    params: control=lambda wildcards, input: '-c ' + input.control if hasattr(input, 'control') else ''
    conda:
        '../Envs/macs.yaml'
    shell:
        'macs2 callpeak -t {input.chip} {params.control} -n IDR/{wildcards.idrprefix} -g {config[genomesize]} --keep-dup all -p 0.01 && cat {output} | sort -k1,1 -k2,2n > {output}.temp && mv {output}.temp {output}'

rule call_broad_idr_peaks:
    input:
        unpack(idr_peak_call_inputs)
    output:
        'IDR/{idrprefix}_peaks.broadPeak'
    params: control=lambda wildcards, input: '-c ' + input.control if hasattr(input, 'control') else ''
    conda:
        '../Envs/macs.yaml'
    shell:
        'macs2 callpeak -t {input.chip} {params.control} -n IDR/{wildcards.idrprefix} -g {config[genomesize]} --keep-dup all -p 0.1 --broad && cat {output} | sort -k1,1 -k2,2n > {output}.temp && mv {output}.temp {output}'

rule pooled_pseudo_peaks:
    input: 
        pool='IDR/{assay}_{tissue}.Pooled_peaks.narrowPeak',
        rep1='IDR/{assay}_{tissue}.Poolr1_peaks.narrowPeak',
        rep2='IDR/{assay}_{tissue}.Poolr2_peaks.narrowPeak'
    output:
        'Peak_Calls/{assay}_{tissue}_IDR.narrowPeak'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        """
        intersectBed -wo -a {input.pool} -b {input.rep1} | 
        awk 'BEGIN{{FS="\\t";OFS="\\t"}}{{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-10 | sort | uniq | 
        intersectBed -wo -a stdin -b {input.rep2} | 
        awk 'BEGIN{{FS="\\t";OFS="\\t"}}{{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-10 | sort | uniq | sort -k1,1 -k2,2n > {output}
        """

rule pooled_pseudo_peaks_broad:
    input: 
        pool='IDR/{assay}_{tissue}.Pooled_peaks.broadPeak',
        rep1='IDR/{assay}_{tissue}.Poolr1_peaks.broadPeak',
        rep2='IDR/{assay}_{tissue}.Poolr2_peaks.broadPeak'
    output:
        'Peak_Calls/{assay}_{tissue}_IDR.broadPeak'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        """
        intersectBed -wo -a {input.pool} -b {input.rep1} | 
        awk 'BEGIN{{FS="\\t";OFS="\\t"}}{{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-9 | sort | uniq | 
        intersectBed -wo -a stdin -b {input.rep2} | 
        awk 'BEGIN{{FS="\\t";OFS="\\t"}}{{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-9 | sort | uniq | sort -k1,1 -k2,2n > {output}
        """
        
rule idr_frip_graph:
    input: 
        bam=lambda wildcards: expand('Aligned_Reads/{library}.bam', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue)), 
        peaks=idr_peak_file
    output: 
        metrics='Metrics/{assay}_{tissue}_IDR_FRiP.txt', 
        figure='Figures/{assay}_{tissue}_IDR_FRiP.png'
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell: 'plotEnrichment --bamfiles {input.bam} --BED {input.peaks} -o {output.figure} --outRawCounts {output.metrics} -e 200 -p {threads}'

##############
# Call Peaks #
##############
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

rule call_atac_peaks:
    input: 
        chip = 'Aligned_Reads/ATAC_{sample}.bam'
    output:
        'Macs2/ATAC_{sample}_peaks.broadPeak' 
    conda:
        '../Envs/macs.yaml'
    threads: 8
    shell:
        'macs2 callpeak -t {input.chip} -f BAM -n Macs2/ATAC_{wildcards.sample} -g {config[genomesize]} -q 0.05 --broad --nomodel --shift -100 --extsize 200'

rule rename_atac_peaks:
    input:
        'Macs2/ATAC_{sample}_peaks.broadPeak'
    output:
        'Peak_Calls/ATAC_{sample}_Peaks.bed'
    shell:
        'cp {input} {output}'

rule fold_enrichment:
    input: 
        'Peak_Calls/{library}_Peaks.bed'
    output: 
        'Macs2/{library}_FoldEnrichment.bdg'
    conda:
        '../Envs/macs.yaml'
    threads: 8
    shell:
        'macs2 bdgcmp -t Macs2/{wildcards.library}_treat_pileup.bdg -c Macs2/{wildcards.library}_control_lambda.bdg -o {output} -m FE -p 0.000001 &&'
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
        'macs2 bdgcmp -t Macs2/{wildcards.library}_treat_pileup.bdg -c Macs2/{wildcards.library}_control_lambda.bdg -o {output} -m logLR -p 0.000001 &&'
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
        bam='Aligned_Reads/{library}.bam',
        bai = 'Aligned_Reads/{library}.bam.bai',
        peaks = 'Peak_Calls/{library}_Peaks.bed'
    output:
        metrics = 'Metrics/{library}_FRiP.txt',
        figure = temp('Metrics/{library}_FRiP.png')
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    script:
        '../Scripts/Calculate_FRiP.py'
    #run:
        #import os
        #if os.path.getsize(input.peaks) > 0:
            #shell('module load python deepTools && plotEnrichment --bamfiles {input.bam} --BED {input.peaks} -o {output.figure} --outRawCounts {output.metrics} -e 200 -p {threads}')
        #else:
            #shell('touch {output.metrics} {output.figure}')
        
rule combined_frip:
    input: 
        bam=lambda wildcards: expand('Aligned_Reads/{library}.bam', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue)), 
        bais = lambda wildcards: expand('Aligned_Reads/{library}.bam.bai', library=libraries(assay=wildcards.assay, tissue=wildcards.tissue)),
        peaks=combined_peak_file
    output: metrics='Metrics/{assay}_{tissue}_FRiP.txt', figure=temp('Metrics/{assay}_{tissue}_FRiP.png')
    conda:
        '../Envs/deeptools.yaml'
    threads: 24
    shell: 'plotEnrichment --bamfiles {input.bam} --BED {input.peaks} -o {output.figure} --outRawCounts {output.metrics} -e 200 -p {threads}'
