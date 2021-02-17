rule Trim_Bedgraph:
    input:
        bdg = 'Macs2/{library}_FoldEnrichment.bdg',
        chromsizes = config['chromsizes']
    output:
        temp('Track_Hub/{library}_FoldEnrichment_Trimmed.bdg')
    conda:
        '../Envs/bdg2bw.yaml'
    shell:
        'bedtools slop -i {input.bdg} -g {input.chromsizes} -b 0 | /home/ckern/bin/bedClip stdin {input.chromsizes} {output}'

rule FoldEnrichment_BigWig:
    input:
        bdg = 'Track_Hub/{library}_FE_Merged.bdg',
        chromsizes = config['chromsizes']
    output:
        'Track_Hub/{library}_FoldEnrichment.bw'
    conda:
        '../Envs/bdg2bw.yaml'
    shell:
        'bedGraphToBigWig {input} {output}'

rule BamCoverage_BigWig:
    input:
        bdg = 'DeepTools/{library}.bw',
    output:
        'Track_Hub/{library}_Coverage.bw'
    shell:
        'ln -sr {input.bdg} {output}'

rule BamCompare_BigWig:
    input:
        bdg = 'DeepTools/{library}_vs_Input.bw',
    output:
        'Track_Hub/{library}_Track.bw'
    shell:
        'ln -sr {input.bdg} {output}'

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
    conda:
        '../Env/bdg2bw.yaml'
    shell:
        'bedToBigBed -type=bed4+1 {input} {output}'

rule BroadPeak_BigBed:
    input:
        peaks = 'Peak_Calls/{assay}_{tissue}_Combined_Peaks.bed_trimmed',
        chromsizes = config['chromsizes']
    output:
        'Track_Hub/{assay}_{tissue}_Broad.bigBed'
    conda:
        '../Env/bdg2bw.yaml'
    shell:
        'bedToBigBed -type=bed4+5 {input} {output}'

rule Temp_Peak_File:
    input:
        peaks = 'Peak_Calls/{library}_Peaks.bed',
        chromsizes = config['chromsizes']
    output:
        temp('Track_Hub/{library}_Temp_Peaks.bed')
    conda:
        '../Envs/bdg2bw.yaml'
    shell:
        'grep -v chrM {input.peaks} | cut -f1,2,3,4 > {input.peaks}.temp &&'
        'bedtools slop -i {input.peaks}.temp -g {input.chromsizes} -b 0 | /home/ckern/bin/bedClip stdin {input.chromsizes} {output} &&'
        'rm {input.peaks}.temp'

rule RepPeak_BigBed:
    input:
        peaks = 'Track_Hub/{library}_Temp_Peaks.bed',
        chromsizes = config['chromsizes']
    output:
        'Track_Hub/{library}_Peaks.bigBed'
    conda:
        '../Env/bdg2bw.yaml'
    shell:
        'bedToBigBed -type=bed4 {input} {output}'

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

rule Link_Bam:
    input:
        'Aligned_Reads/{library}.bam'
    output:
        'Track_Hub/{library}.bam'
    shell:
        'ln -sr {input} {output} && ln -sr {input}.bai {output}.bai'

rule Create_TrackDB:
    input:
        #fes = expand('Track_Hub/{library}_Track.bw', library=peak_libraries())
        #fes = expand('Track_Hub/{library}_FoldEnrichment.bw', library=peak_libraries()),
        fes = expand('Track_Hub/{library}_Coverage.bw', library=libraries(skip_input=False)),
        #peaks = expand('Track_Hub/{library}_Peaks.bigBed', library=libraries(peak_assay_only=True))
        #fes = expand('Track_Hub/{library}_Pileup.bw', library=peak_libraries()),
        #narrow = expand('Track_Hub/{library}_Combined_Peaks.bigBed', library=combined_libraries(config['narrow_peaks'])),
        #broad = expand('Track_Hub/{library}_Combined_Peaks.bigBed', library=combined_libraries(config['broad_peaks']))
        #segs = expand('Track_Hub/{tissue}_Chromatin_State_{size}.bigBed', tissue=tissues_for_assay(assay='H3K4me3'), size=['Dense'])
        #bams = expand('Track_Hub/{library}.bam', library=libraries(skip_input=False))
    output:
        trackdb = 'Track_Hub/trackDb.txt'
    script:
        '../Scripts/Make_TrackDB.py'
