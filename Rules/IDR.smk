
def idr_peak_call_inputs(wildcards):
    prefix = config['tempdir'] + '/{}_{}_IDR/'.format(wildcards.assay, wildcards.tissue)
    inputs = {'chip': prefix + '{}.tagAlign.gz'.format(wildcards.idrprefix)}
    if assay not in config['no_input']:
        #if wildcards.idrprefix.startswith('Rep1'):
            #inputs['control'] = prefix + 'Inp1.tagAlign.gz'
        #elif wildcards.idrprefix.startswith('Rep2'):
            #inputs['control'] = prefix + 'Inp2.tagAlign.gz'
        #elif wildcards.idrprefix.startswith('RepPool'):
        inputs['control'] = prefix + 'InpPool.tagAlign.gz'
    return inputs

def copy_inputs(wildcards):
    inputs = {}
    for i, rep in enumerate(libraries(assay=wildcards.assay, tissue=wildcards.tissue)):
        inputs['rep' + str(i + 1) + 'chip'] = 'Aligned_Reads/{}.tagAlign.gz'.format(rep)
        inputs['rep' + str(i + 1) + 'input'] = control_library(rep, format='tagAlign')
    return inputs

rule copy_alignments:
    input:
        unpack(copy_inputs)
    output:
        rep1chip = temp(config['tempdir'] + '/{assay}_{tissue}_IDR/Rep1.tagAlign.gz'),
        rep2chip = temp(config['tempdir'] + '/{assay}_{tissue}_IDR/Rep2.tagAlign.gz'),
        rep1input = temp(config['tempdir'] + '/{assay}_{tissue}_IDR/Inp1.tagAlign.gz'),
        rep2input = temp(config['tempdir'] + '/{assay}_{tissue}_IDR/Inp2.tagAlign.gz')
    #group: 'idr'
    shell:
        """
        cp {input.rep1chip} {output.rep1chip}
        cp {input.rep2chip} {output.rep2chip}
        cp {input.rep1input} {output.rep1input}
        cp {input.rep2input} {output.rep2input}
        """
        

#Splits a set of alignments into random halves
rule create_pseudoreplicates:
    input: 
        config['tempdir'] + '/{assay}_{tissue}_IDR/{idrprefix}.tagAlign.gz'
    output: 
        pr1 = temp(config['tempdir'] + '/{assay}_{tissue}_IDR/{idrprefix}Pr1.tagAlign.gz'), 
        pr2 = temp(config['tempdir'] + '/{assay}_{tissue}_IDR/{idrprefix}Pr2.tagAlign.gz')
    params:
        prefix = config['tempdir'] + '/{assay}_{tissue}_IDR'
    #group: 'idr'
    shell:
        """
        nlines=$( zcat {input} | wc -l )
        nlines=$(( (nlines + 1) / 2 ))
        zcat {input} | shuf --random-source={input} | split -d -l ${{nlines}} - {params.prefix}/{wildcards.idrprefix}
        gzip -nc "{params.prefix}/{wildcards.idrprefix}00" > {output.pr1}
        rm "{params.prefix}/{wildcards.idrprefix}00"
        gzip -nc "{params.prefix}/{wildcards.idrprefix}01" > {output.pr2}
        rm "{params.prefix}/{wildcards.idrprefix}01"
        """

#Combines sets of alignments into a single pool
rule create_pooled:
    input: 
        config['tempdir'] + '/{assay}_{tissue}_IDR/{type}1.tagAlign.gz',
        config['tempdir'] + '/{assay}_{tissue}_IDR/{type}2.tagAlign.gz'
    output: 
        temp(config['tempdir'] + '/{assay}_{tissue}_IDR/{type}Pool.tagAlign.gz')
    #group: 'idr'
    shell:
        'zcat {input} | gzip -nc > {output}'

#Call peaks using macs2
rule call_idr_peaks:
    input:
        unpack(idr_peak_call_inputs)
    output:
        temp(config['tempdir'] + '/{assay}_{tissue}_IDR/{idrprefix}_peaks.narrowPeak')
    params: 
        control = lambda wildcards, input: '-c ' + input.control if hasattr(input, 'control') else ''
    #group: 'idr'
    conda:
        '../Envs/macs.yaml'
    shell:
        'macs2 callpeak -t {input.chip} {params.control} -n {config[tempdir]}/{wildcards.assay}_{wildcards.tissue}_IDR/{wildcards.idrprefix} -g {config[genomesize]} --keep-dup all -q 0.05 && cat {output} | sort -k1,1 -k2,2n > {output}.temp && mv {output}.temp {output}'

#Call peaks using spp
#rule call_idr_peaks:
    #input:
        #unpack(idr_peak_call_inputs)
    #output:
        #temp(config['tempdir'] + '/{assay}_{tissue}_IDR/{idrprefix}_peaks.narrowPeak')
    #params: 
        #control = lambda wildcards, input: '-i=' + input.control if hasattr(input, 'control') else ''
    #group: 'idr'
    #threads: 12
    #conda:
        #'../Envs/spp.yaml'
    #shell:
        #'run_spp.R -c={input.chip} {params.control} -npeak=300000 -odir={config[tempdir]}/{assay}_{tissue}_IDR/ -rf -savr -p=12'


def idr_inputs(wildcards):
    prefix = config['tempdir'] + '/{}_{}_IDR/'.format(wildcards.assay, wildcards.tissue)
    if wildcards.type == 'TruePooled':
        return {
            'pool': prefix + 'RepPool_peaks.narrowPeak',
            'rep1': prefix + 'Rep1_peaks.narrowPeak',
            'rep2': prefix + 'Rep2_peaks.narrowPeak'
        }
    elif wildcards.type == 'PseudoPooled':
        return {
            'pool': prefix + 'RepPool_peaks.narrowPeak',
            'rep1': prefix + 'RepPoolPr1_peaks.narrowPeak',
            'rep2': prefix + 'RepPoolPr2_peaks.narrowPeak'
        }
    elif wildcards.type == 'Rep1SelfPr':
        return {
            'pool': prefix + 'Rep1_peaks.narrowPeak',
            'rep1': prefix + 'Rep1Pr1_peaks.narrowPeak',
            'rep2': prefix + 'Rep1Pr2_peaks.narrowPeak'
        }
    elif wildcards.type == 'Rep2SelfPr':
        return {
            'pool': prefix + 'Rep2_peaks.narrowPeak',
            'rep1': prefix + 'Rep2Pr1_peaks.narrowPeak',
            'rep2': prefix + 'Rep2Pr2_peaks.narrowPeak'
        }
    else:
        raise TypeError
                
rule idr_threshold_peaks:
    input: 
        unpack(idr_inputs)
    output:
        temp(config['tempdir'] + '/{assay}_{tissue}_IDR/{type}.bed')
    #group: 'idr'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        """
        intersectBed -wo -a {input.pool} -b {input.rep1} | 
        awk 'BEGIN{{FS="\\t";OFS="\\t"}}{{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-10 | sort | uniq | 
        intersectBed -wo -a stdin -b {input.rep2} | 
        awk 'BEGIN{{FS="\\t";OFS="\\t"}}{{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-10 | sort | uniq | sort -k1,1 -k2,2n > {output}
        """

#rule run_idr:
    #input: 
        #unpack(idr_inputs)
    #output:
        #temp(config['tempdir'] + '/{assay}_{tissue}_IDR/{type}_unfiltered.bed')
    #group: 'idr'
    #conda:
        #'../Envs/idr.yaml'
    #shell:
        #'idr --samples {input.rep1} {input.rep2} --peak-list {input.pool} --input-file-type narrowPeak --output-file {output} --rank signal.value --soft-idr-threshold 0.05 --use-best-multisummit-IDR'

#rule filter_idr_peaks:
    #input:
        #rules.run_idr.output
    #output:
        #'Peak_Calls/{assay}_{tissue}_IDR_{type}.bed'
    #group: 'idr'
    #shell:
        #"""awk 'BEGIN{{OFS="\t"}} $12>=1.301 {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {input} | sort | uniq | sort -k7n,7n > {output}"""

rule get_stable_peaks:
    input:
        config['tempdir'] + '/{assay}_{tissue}_IDR/TruePooled.bed',
        config['tempdir'] + '/{assay}_{tissue}_IDR/PseudoPooled.bed'
    output:
        stable = 'Peak_Calls/{assay}_{tissue}_Stable_Peaks.bed', 
        cons = 'Peak_Calls/{assay}_{tissue}_Conservative_Peaks.bed', 
        pooled = 'Peak_Calls/{assay}_{tissue}_TruePooled_Peaks.bed'
    #group: 'idr'
    conda:
        '../Envs/bedtools.yaml'
    shell:
        'cat {input} | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,5,6,7,8,9,10 -o distinct,mean,distinct,mean,mean,mean,mean > {output.stable} && '
        'bedtools intersect -a {input[0]} -b {input[1]} -u > {output.cons} && '
        'cp {input[0]} {output.pooled}'

rule idr_report:
    input:
        pooled_peaks = config['tempdir'] + '/{assay}_{tissue}_IDR/RepPool_peaks.narrowPeak',
        true_pooled = config['tempdir'] + '/{assay}_{tissue}_IDR/TruePooled.bed',
        pseudo_pooled = config['tempdir'] + '/{assay}_{tissue}_IDR/PseudoPooled.bed',
        stable_peaks = 'Peak_Calls/{assay}_{tissue}_Stable_Peaks.bed',
        rep1_self = config['tempdir'] + '/{assay}_{tissue}_IDR/Rep1SelfPr.bed',
        rep2_self = config['tempdir'] + '/{assay}_{tissue}_IDR/Rep2SelfPr.bed'
    output:
        json = 'Metrics/{assay}_{tissue}_IDR_Stats.json'
    #group: 'idr'
    threads: 24 #This is just to get more memory REMOVE LATER
    script:
        '../Scripts/Get_IDR_Stats.py'

def idr_summary_inputs(wildcards):
    for tissue in tissues_for_assay(assay=wildcards.assay):
        yield 'Metrics/{}_{}_IDR_Stats.json'.format(wildcards.assay, tissue)

rule idr_summary:
    input:
        #expand('Metrics/{combined_library}_IDR_Stats.json', combined_library=combined_libraries(assay=assay))
        idr_summary_inputs
    output:
        'Tables/{assay}_IDR_Summary.txt'
    shell:
        'touch {output}'
