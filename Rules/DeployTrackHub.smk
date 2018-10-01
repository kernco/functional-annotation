rule FoldEnrichment_BigWig:
    input:
        bdg = 'Macs2/{library}_FoldEnrichment.bdg',
        chromsizes = config['chromsizes']
    output:
        'Track_Hub/{library}_FoldEnrichment.bw'
    shell:
        'bedGraphToBigWig {input} {output}'

rule BamCoverage_BigWig:
    input:
        bdg = 'DeepTools/{library}.bw',
        chromsizes = config['chromsizes']
    output:
        'Track_Hub/{library}_Coverage.bw'
    shell:
        'ln {input.bdg} {output}'

rule Sort_TreatPileup:
    input:
        bdg = 'Macs2/{library}_treat_pileup.bdg'
    output:
        temp('Track_Hub/{library}_Sorted.bdg')
    shell:
        'sort -k1,1 -k2,2n {input} > {output}'

rule TreatPileup_BigWig:
    input:
        bdg = 'Track_Hub/{library}_Sorted.bdg_trimmed',
        chromsizes = config['chromsizes']
    output:
        'Track_Hub/{library}_Pileup.bw'
    shell:
        'bedGraphToBigWig {input} {output}'
        
rule Trim_Bed:
    input:
        bed = '{prefix}',
        chromsizes = config['chromsizes']
    output:
        trimmed = temp('{prefix}_trimmed')
    script:
        '../Scripts/TrimBedToChromosomes.py'

rule Rescore_NarrowPeak:
    input:
        narrowpeak = 'Peak_Calls/{prefix}.bed'
    output:
        scaledbed = 'Track_Hub/{prefix}.fixed.narrowPeak'
    script:
        '../Scripts/Rescore_Peaks.py'
        
rule NarrowPeak_BigBed:
    input:
        peaks = 'Track_Hub/{assay}_{tissue}_Combined_Peaks.fixed.narrowPeak',
        chromsizes = config['chromsizes']
    output:
        'Track_Hub/{assay}_{tissue}_Combined_Peaks.bigBed'
    shell:
        'bedToBigBed -type=bed4+1 {input} {output}'

rule BroadPeak_BigBed:
    input:
        peaks = 'Peak_Calls/{assay}_{tissue}_Combined_Peaks.bed_trimmed',
        chromsizes = config['chromsizes']
    output:
        'Track_Hub/{assay}_{tissue}_Broad.bigBed'
    shell:
        'bedToBigBed -type=bed4+5 {input} {output}'

rule Segmentation_Dense_BigBed:
    input:
        'ChromHMM/Model_Bam_15/{tissue}_15_dense.bed'
    output:
        'Track_Hub/{tissue}_Chromatin_State_Dense.bigBed'
    shell:
        'tail -n +2 {input} | sort -k1,1 -k2,2n > {input}_temp && bedToBigBed -type=bed6+3 {input}_temp {config[chromsizes]} {output} && rm {input}_temp'

rule Segmentation_Expanded_BigBed:
    input:
        'ChromHMM/Model_Bam_15/{tissue}_15_expanded.bed'
    output:
        'Track_Hub/{tissue}_Chromatin_State_Expanded.bigBed'
    shell:
        'tail -n +3 {input} | sort -k1,1 -k2,2n > {input}_temp && bedToBigBed {input}_temp {config[chromsizes]} {output} && rm {input}_temp'

rule Create_TrackDB:
    input:
        #fes = expand('Track_Hub/{library}_FoldEnrichment.bw', library=peak_libraries()),
        #fes = expand('Track_Hub/{library}_Coverage.bw', library=peak_libraries()),
        fes = expand('Track_Hub/{library}_Pileup.bw', library=peak_libraries()),
        narrow = expand('Track_Hub/{library}_Combined_Peaks.bigBed', library=combined_libraries(config['narrow_peaks'])),
        broad = expand('Track_Hub/{library}_Combined_Peaks.bigBed', library=combined_libraries(config['broad_peaks']))
        #segs = expand('Track_Hub/{tissue}_Chromatin_State_{size}.bigBed', tissue=tissues_for_assay(assay='H3K4me3'), size=['Dense'])
    output:
        trackdb = 'Track_Hub/trackDb.txt'
    script:
        '../Scripts/Make_TrackDB.py'
