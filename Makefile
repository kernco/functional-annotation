
.SECONDARY:
.DELETE_ON_ERROR:
.SUFFIXES:

DIRS = Trimmed_Reads Bwa_Output Tophat_Output Aligned_Reads Peak_Calls Results DeepTools Cufflinks_Output Bed_Files Temp_Files Chromatin_Model LncRNA
PIPELINE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

###############################
# Load Pipeline Configuration #
###############################
define Load_Setting
$1 := $(shell grep $1 pipeline-config.txt | cut -d'=' -f2)
endef
SETTINGS := ASSEMBLY ANNOTATION BROAD_PEAKS_WITH_INPUT NARROW_PEAKS_WITH_INPUT BROAD_PEAKS_NO_INPUT GENOME_SIZE FRAGMENT_LENGTH
$(foreach setting,$(SETTINGS),$(eval $(call Load_Setting,$(setting))))

ASSEMBLY_BASE := $(basename $(ASSEMBLY))
ANNOTATION_BASE := $(basename $(ANNOTATION))
PEAK_ASSAYS := $(BROAD_PEAKS_NO_INPUT) $(NARROW_PEAKS_WITH_INPUT) $(BROAD_PEAKS_WITH_INPUT)
BWA_ASSAYS := $(PEAK_ASSAYS) Input
ASSAYS_WITH_INPUT := $(NARROW_PEAKS_WITH_INPUT) $(BROAD_PEAKS_WITH_INPUT)
BROAD_PEAKS := $(BROAD_PEAKS_NO_INPUT) $(BROAD_PEAKS_WITH_INPUT)
ASSAYS := RNASeq $(BWA_ASSAYS)

#############################################
# Get Tissues and Replicates for Each Assay #
#############################################
define Get_Libraries
$(eval temp := $(patsubst Raw_Reads/%.fq.gz,%,$(wildcard Raw_Reads/$1_*.fq.gz)))
$(eval temp := $(subst _R1,,$(temp)))
$(eval temp := $(subst _R2,,$(temp)))
$(eval $1_LIBRARIES :=)
$(foreach lib,$(temp),$(if $(filter $(lib),$($1_LIBRARIES)),,$(eval $1_LIBRARIES += $(lib))))
endef
$(foreach assay,$(ASSAYS),$(eval $(call Get_Libraries,$(assay))))

#Duplicates will be in here, but $^ omits duplicate prerequisites
define Add_Samples
$(foreach sample,$($1_LIBRARIES:$1_%=%),$(if $(filter $(sample),$(SAMPLES)),,$(eval SAMPLES += $(sample))))
endef
$(foreach assay,$(ASSAYS),$(eval $(call Add_Samples,$(assay))))

define Get_Tissues
$(foreach lib,$($1_LIBRARIES:$1_%=%),$(if $(filter $(word 1,$(subst _, ,$(lib))),$($1_TISSUES)),,$(eval $1_TISSUES += $(word 1,$(subst _, ,$(lib))))))
endef
$(foreach assay,$(ASSAYS),$(eval $(call Get_Tissues,$(assay))))

define Get_All_Tissues
$(foreach tissue,$($1_TISSUES),$(if $(filter $(tissue),$(TISSUES)),,$(eval TISSUES += $(tissue))))
endef
$(foreach assay,$(ASSAYS),$(eval $(call Get_All_Tissues,$(assay))))

define Get_Replicates
$(foreach lib,$(filter $2%,$($1_LIBRARIES:$1_%=%)),$(if $(filter $(word 2,$(subst _, ,$(lib))),$($1_$2_REPLICATES)),,$(eval $1_$2_REPLICATES += $(word 2,$(subst _, ,$(lib))))))
endef
$(foreach assay,$(ASSAYS),$(foreach tissue,$($(assay)_TISSUES),$(eval $(call Get_Replicates,$(assay),$(tissue)))))

define Get_All_Replicates
$(foreach rep,$($1_$2_REPLICATES),$(if $(filter $(rep),$($1_REPLICATES)),,$(eval $1_REPLICATES += $(rep))))
endef
$(foreach assay,$(ASSAYS),$(foreach tissue,$($(assay)_TISSUES),$(eval $(call Get_All_Replicates,$(assay),$(tissue)))))

define Get_Assays
$(eval temp := $(patsubst Raw_Reads/%_$1.fq.gz,%,$(subst _R2,,$(subst _R1,,$(wildcard Raw_Reads/*_$1*.fq.gz)))))
$(eval $1_ASSAYS :=)
$(foreach assay,$(temp),$(if $(filter $(assay),$($1_ASSAYS)),,$(eval $1_ASSAYS += $(assay))))
endef
$(foreach sample,$(SAMPLES),$(eval $(call Get_Assays,$(sample))))

######################
# High-Level Targets #
######################
PHONIES := all AlignAll
all : AlignAll CallAllPeaks

CallAllPeaks : $(foreach assay,$(PEAK_ASSAYS),Call$(assay)Peaks)

define Peak_Template
ifdef $1_LIBRARIES
Call$1Peaks : Results/$1_Peak_Heatmap.png Results/$1_Peak_Summary.txt
else
Call$1Peaks :
endif
endef
$(foreach assay,$(NARROW_PEAKS_WITH_INPUT),$(eval $(call Peak_Template,$(assay),narrow)))
$(foreach assay,$(BROAD_PEAKS),$(eval $(call Peak_Template,$(assay),broad)))


AlignAll : $(foreach assay,$(ASSAYS),Align$(assay)) $(foreach sample,$(SAMPLES),Results/$(sample)_Fingerprint.png Results/$(sample)_Coverage.png) ProfileGraphs 

ifdef RNASeq_LIBRARIES
ProfileGraphs : $(foreach assay,$(PEAK_ASSAYS),$(foreach lib,$($(assay)_LIBRARIES),Results/$(lib)_Profile.png Results/$(lib)_Heatmap.png))
else
ProfileGraphs :
endif

define Align_Targets
ifdef $1_LIBRARIES
Align$1 : Results/$1_Alignment_Summary.txt Results/$1_Pearson_Correlation.png Results/$1_Spearman_Correlation.png $(if $(filter $1,$(ASSAYS_WITH_INPUT)),Results/$1_vs_Input_Pearson_Correlation.png,) $(if $(filter $1,$(ASSAYS_WITH_INPUT)),Results/$1_vs_Input_Spearman_Correlation.png,) $(if $(filter $1,$(BWA_ASSAYS)),Results/$1_Quality_Metrics.txt)
else
Align$1 :
endif
endef
$(foreach assay,$(ASSAYS),$(eval $(call Align_Targets,$(assay))))
PHONIES += $(foreach assay,$(ASSAYS),Align$(assay))

############
# ChromHMM #
############
Chromatin_Model/Model_% : Chromatin_Model/Binarized_Data
	java -mx10000M -jar ~/ChromHMM/ChromHMM.jar LearnModel -p $(THREADS) -r 5000 -l Chromatin_Model/Chromosome_Lengths.txt $< $@ $* galGal_6tissues

Chromatin_Model/Binarized_Data : Chromatin_Model/Chromosome_Lengths.txt Chromatin_Model/Tissue_Marks.txt $(foreach assay,$(PEAK_ASSAYS) RNASeq Input,Chromatin_Model/$(assay).merged.bam)
	java -mx10000M -jar ~/ChromHMM/ChromHMM.jar BinarizeBam Chromatin_Model/Chromosome_Lengths.txt Chromatin_Model Chromatin_Model/Tissue_Marks.txt $@

Chromatin_Model/Tissue_Marks.txt : | Chromatin_Model
	rm -f $@
	$(foreach assay,$(ASSAYS_WITH_INPUT),echo "Combined\t$(assay)\t$(assay).merged.bam\tInput.merged.bam" >> $@;)
	$(foreach assay,$(BROAD_PEAKS_NO_INPUT) RNASeq,echo "Combined\t$(assay)\t$(assay).merged.bam" >> $@;)

#$(foreach assay,$(ASSAYS_WITH_INPUT),$(foreach lib,$(subst $(assay)_,,$($(assay)_LIBRARIES)),echo "$(lib)\t$(assay)\t$(assay)_$(lib).bam\tInput_$(lib).bam" >> $@;))
#$(foreach assay,$(BROAD_PEAKS_NO_INPUT),$(foreach lib,$(subst $(assay)_,,$($(assay)_LIBRARIES)),echo "$(lib)\t$(assay)\t$(assay)_$(lib).bam" >> $@;))
	
Chromatin_Model/Chromosome_Lengths.txt : | Chromatin_Model
	python $(PIPELINE_DIR)/Get_Chromosome_Lengths.py $(ASSEMBLY) $(TISSUES) > $@

define Merge_Bams
Chromatin_Model/$1.merged.bam : $(foreach lib,$($1_LIBRARIES),Chromatin_Model/$(lib).bam)
	samtools merge $$@ $$^
endef
$(foreach assay,$(ASSAYS),$(eval $(call Merge_Bams,$(assay))))
#$(foreach assay,$(BROAD_PEAKS_NO_INPUT),$(eval $(call Merge_Bams,$(assay))))

Chromatin_Model/%.bam : Aligned_Reads/%.bam | Chromatin_Model
	samtools view -h $< | python $(PIPELINE_DIR)/Prep_Bam_For_ChromHMM.py $* | samtools view -b - > $@


######################
# Fingerprint Graphs #
######################
define Fingerprint_Template
$(eval deps := $(foreach assay,$(BWA_ASSAYS),$(wildcard Raw_Reads/$(assay)_$1.fq.gz)))
$(eval deps += $(foreach assay,$(BWA_ASSAYS),$(wildcard Raw_Reads/$(assay)_$1_R1.fq.gz)))
$(eval deps := $(subst _R1,,$(deps)))
$(eval deps := $(patsubst Raw_Reads/%.fq.gz,Aligned_Reads/%.bam.bai,$(deps)))
Results/$1_Fingerprint.png : $(deps) | Results $(foreach assay,$(BWA_ASSAYS),Align$(assay))
	plotFingerprint -b $$(^:.bai=) -plot $$@ --labels $$(^:Aligned_Reads/%_$1.bam.bai=%) --plotTitle="$(subst _, ,$1)" -p=$(THREADS) --extendReads=$(FRAGMENT_LENGTH) --skipZeros
endef
$(foreach sample,$(SAMPLES),$(eval $(call Fingerprint_Template,$(sample))))

##################
# Profile Graphs #
##################
ifdef RNASeq_LIBRARIES
define Plot_Profile
Results/$1_Profile.png : DeepTools/$1_GeneExp.mat.gz | Results
	plotProfile -m $$< -out $$@ --plotTitle "$(subst _, ,$1)"

Results/$1_Heatmap.png : DeepTools/$1_GeneExp.mat.gz | Results
	plotHeatmap -m $$< -out $$@ --plotTitle "$(subst _, ,$1)"
endef
$(foreach assay,$(PEAK_ASSAYS),$(foreach lib,$($(assay)_LIBRARIES),$(eval $(call Plot_Profile,$(lib)))))

define Profile_Matrix_NoInput
DeepTools/$1_GeneExp.mat.gz : DeepTools/$1.bw Bed_Files/RNASeq_$2/Expressed_Genes.bed Bed_Files/RNASeq_$2/Repressed_Genes.bed | DeepTools
	computeMatrix reference-point --referencePoint TSS -S DeepTools/$1.bw -R Bed_Files/RNASeq_$2/Expressed_Genes.bed Bed_Files/RNASeq_$2/Repressed_Genes.bed -b 3000 -a 10000 --skipZeros -out $$@ -p $(THREADS)
endef
$(foreach assay,$(BROAD_PEAKS_NO_INPUT),$(foreach lib,$($(assay)_LIBRARIES),$(eval $(call Profile_Matrix_NoInput,$(lib),$(subst $(assay)_,,$(lib))))))

define Profile_Matrix_Input
DeepTools/$1_GeneExp.mat.gz : DeepTools/$1_vs_Input.bw Bed_Files/RNASeq_$2/Expressed_Genes.bed Bed_Files/RNASeq_$2/Repressed_Genes.bed | DeepTools
	computeMatrix reference-point --referencePoint TSS -S DeepTools/$1_vs_Input.bw -R Bed_Files/RNASeq_$2/Expressed_Genes.bed Bed_Files/RNASeq_$2/Repressed_Genes.bed -b 3000 -a 10000 --skipZeros -out $$@ -p $(THREADS)
endef
$(foreach assay,$(ASSAYS_WITH_INPUT),$(foreach lib,$($(assay)_LIBRARIES),$(eval $(call Profile_Matrix_Input,$(lib),$(subst $(assay)_,,$(lib))))))
endif

################################
# Alignment Correlation Graphs #
################################
define Correlation_Graph
Results/$1_$2_Correlation.png : DeepTools/$1_MultiBigwigSummary.npz | Results
	plotCorrelation -in $$< --corMethod $3 --skipZeros --plotTitle "$2 Correlation of $1 Read Counts" --whatToPlot heatmap --plotNumbers -o $$@ --removeOutliers
endef
$(foreach assay,$(ASSAYS),$(eval $(call Correlation_Graph,$(assay),Pearson,pearson)))
$(foreach assay,$(ASSAYS),$(eval $(call Correlation_Graph,$(assay),Spearman,spearman)))
$(foreach assay,$(ASSAYS_WITH_INPUT),$(eval $(call Correlation_Graph,$(assay)_vs_Input,Pearson,pearson)))
$(foreach assay,$(ASSAYS_WITH_INPUT),$(eval $(call Correlation_Graph,$(assay)_vs_Input,Spearman,spearman)))

ifdef RNASeq_LIBRARIES
DeepTools/RNASeq_MultiBigwigSummary.npz : $(foreach lib,$(RNASeq_LIBRARIES),DeepTools/$(lib).bw) | DeepTools
	multiBigwigSummary bins -p=$(THREADS) -b $^ -out $@ --labels $(subst RNASeq_,,$(RNASeq_LIBRARIES)) 

DeepTools/RNASeq_%.bw : Aligned_Reads/RNASeq_%.bam Aligned_Reads/RNASeq_%.bam.bai | DeepTools
	bamCoverage -p=$(THREADS) -b $< -o $@ -of=bigwig --normalizeUsingRPKM
endif

define MinusInput_MultiBigwig_Template
DeepTools/$1_vs_Input_MultiBigwigSummary.npz : $(foreach lib,$($1_LIBRARIES),DeepTools/$(lib)_vs_Input.bw) | DeepTools
	multiBigwigSummary bins -p=$(THREADS) -b $$^ -out $$@ --labels $(subst $1_,,$($1_LIBRARIES)) --binSize=1000
endef
$(foreach assay,$(ASSAYS_WITH_INPUT),$(eval $(call MinusInput_MultiBigwig_Template,$(assay))))

define MultiBigwig_Rule
DeepTools/$1_MultiBigwigSummary.npz : $(foreach lib,$($1_LIBRARIES),DeepTools/$(lib).bw) | DeepTools
	multiBigwigSummary bins -p=$(THREADS) -b $$^ -out $$@ --labels $(subst $1_,,$($1_LIBRARIES)) --binSize=1000
endef
$(foreach assay,$(BWA_ASSAYS),$(eval $(call MultiBigwig_Rule,$(assay))))

define Vs_Input_Bigwig
DeepTools/$1_$2_vs_Input.bw : Aligned_Reads/$1_$2.bam Aligned_Reads/Input_$2.bam Aligned_Reads/$1_$2.bam.bai Aligned_Reads/Input_$2.bam.bai
	bamCompare -b1 $$< -b2 $$(word 2,$$^) -o $$@ -p=$(THREADS) --ignoreDuplicates --extendReads=$(FRAGMENT_LENGTH)
endef
$(foreach assay,$(ASSAYS_WITH_INPUT),$(foreach sample,$($(assay)_LIBRARIES:$(assay)_%=%),$(eval $(call Vs_Input_Bigwig,$(assay),$(sample)))))

DeepTools/%.bw : Aligned_Reads/%.bam Aligned_Reads/%.bam.bai | DeepTools
	bamCoverage -p=$(THREADS) -b $< -o $@ -of=bigwig --normalizeUsingRPKM --ignoreDuplicates --extendReads=$(FRAGMENT_LENGTH)

#############################
# Alignment Coverage Graphs #
#############################
define Coverage_Template
$(eval deps := $(foreach assay,$(ASSAYS),$(wildcard Raw_Reads/$(assay)_$1.fq.gz)))
$(eval deps += $(foreach assay,$(ASSAYS),$(wildcard Raw_Reads/$(assay)_$1_R1.fq.gz)))
$(eval deps := $(subst _R1,,$(deps)))
$(eval deps := $(patsubst Raw_Reads/%.fq.gz,Aligned_Reads/%.bam.bai,$(deps)))
Results/$1_Coverage.png : $(deps) | Results 
	plotCoverage -b $$(^:.bai=) -o $$@ --labels $$(^:Aligned_Reads/%_$1.bam.bai=%) --plotTitle="$(subst _, ,$1) Library Coverage" -p=$(THREADS) --extendReads=$(FRAGMENT_LENGTH) --ignoreDuplicates
endef
$(foreach sample,$(SAMPLES),$(eval $(call Coverage_Template,$(sample))))

############################
# Alignment Summary Tables #
############################
define Align_Summary
Results/$1_Alignment_Summary.txt : $(foreach lib,$($1_LIBRARIES),Aligned_Reads/$(lib).bam) | Results
	python $(PIPELINE_DIR)/Make_Alignment_Summary.py $($1_LIBRARIES) > $$@
endef
$(foreach assay,$(ASSAYS),$(eval $(call Align_Summary,$(assay))))

define Quality_Metrics
Results/$1_Quality_Metrics.txt : $(foreach lib,$($1_LIBRARIES),Aligned_Reads/$(lib).bam) | Results
	python $(PIPELINE_DIR)/Library_Quality_Metrics.py "$($1_LIBRARIES)" $(THREADS) > $$@
endef
$(foreach assay,$(BWA_ASSAYS),$(eval $(call Quality_Metrics,$(assay))))

######################
# Combine Replicates #
######################
define Overlap_Peaks
ifneq ($(strip $(word 2,$($3_REPLICATES))),)
Peak_Calls/$3_Combined_Peaks.bed : | Peak_Calls/$1 Peak_Calls/$2
	bedtools intersect -a Peak_Calls/$1/macs2_peaks.bed -b Peak_Calls/$2/macs2_peaks.bed -u -sorted -f 0.5 > Peak_Calls/$3_temp
	bedtools intersect -a Peak_Calls/$2/macs2_peaks.bed -b Peak_Calls/$1/macs2_peaks.bed -u -sorted -f 0.5 >> Peak_Calls/$3_temp
	sort -k1,1 -k2,2n -o Peak_Calls/$3_temp Peak_Calls/$3_temp
	bedtools merge -i Peak_Calls/$3_temp > $$@
	rm Peak_Calls/$3_temp
else
Peak_Calls/$3_Combined_Peaks.bed : | Peak_Calls/$1
	ln Peak_Calls/$1/macs2_peaks.bed $$@
endif
endef
$(foreach assay,$(BROAD_PEAKS),$(foreach tissue,$($(assay)_TISSUES),$(eval $(call Overlap_Peaks,$(assay)_$(tissue)_$(word 1,$($(assay)_$(tissue)_REPLICATES)),$(assay)_$(tissue)_$(word 2,$($(assay)_$(tissue)_REPLICATES)),$(assay)_$(tissue)))))

define Combine_Peaks
ifneq ($(strip $(word 2,$($3_REPLICATES))),)
Peak_Calls/$3_Combined_Peaks.bed : Peak_Calls/$1/Peaks_Validated_by_Replicate.bed Peak_Calls/$2/Peaks_Validated_by_Replicate.bed
	cat $$^ | sort -k1,1 -k2,2n > Peak_Calls/$3_temp
	if [ -s Peak_Calls/$3_temp ]; then bedtools merge -i Peak_Calls/$3_temp > $$@; else touch $$@; fi
	rm Peak_Calls/$3_temp
else
Peak_Calls/$3_Combined_Peaks.bed : | Peak_Calls/$1
	ln Peak_Calls/$1/macs2_peaks.bed $$@
endif
endef
$(foreach assay,$(NARROW_PEAKS_WITH_INPUT),$(foreach tissue,$($(assay)_TISSUES),$(eval $(call Combine_Peaks,$(assay)_$(tissue)_$(word 1,$($(assay)_$(tissue)_REPLICATES)),$(assay)_$(tissue)_$(word 2,$($(assay)_$(tissue)_REPLICATES)),$(assay)_$(tissue)))))

Peak_Calls/%/Peaks_Validated_by_Replicate.bed : Peak_Calls/%/Summits_Validated_by_Replicate.bed
	python $(PIPELINE_DIR)/Summits2Peaks.py $< Peak_Calls/$*/macs2_peaks.bed > $@

define Validate_Summits
ifneq ($(strip $(word 3,$(subst _, ,$2))),)
Peak_Calls/$1/Summits_Validated_by_Replicate.bed : Peak_Calls/$2/FoldEnrichment.bdg Peak_Calls/$2/LogLikelihoodRatio.bdg | Peak_Calls/$1
	bedtools intersect -a Peak_Calls/$1/macs2_summits.bed -b $$< -wa -wb -sorted | bedtools intersect -a stdin -b $$(word 2,$$^) -wa -wb -sorted | cut -f1,2,3,4,9,13 | awk '$$$$5 > 1 && $$$$6 > 0' > $$@
endif
endef
$(foreach assay,$(NARROW_PEAKS_WITH_INPUT),$(foreach tissue,$($(assay)_TISSUES),$(eval $(call Validate_Summits,$(assay)_$(tissue)_$(word 1,$($(assay)_$(tissue)_REPLICATES)),$(assay)_$(tissue)_$(word 2,$($(assay)_$(tissue)_REPLICATES))))))
$(foreach assay,$(NARROW_PEAKS_WITH_INPUT),$(foreach tissue,$($(assay)_TISSUES),$(eval $(call Validate_Summits,$(assay)_$(tissue)_$(word 2,$($(assay)_$(tissue)_REPLICATES)),$(assay)_$(tissue)_$(word 1,$($(assay)_$(tissue)_REPLICATES))))))

Peak_Calls/%/FoldEnrichment.bdg : | Peak_Calls/%
	macs2 bdgcmp -t Peak_Calls/$*/macs2_treat_pileup.bdg -c Peak_Calls/$*/macs2_control_lambda.bdg -o $@ -m FE -p 0.00001
	sort -k1,1 -k2,2n -o $@ $@

Peak_Calls/%/LogLikelihoodRatio.bdg : | Peak_Calls/%
	macs2 bdgcmp -t Peak_Calls/$*/macs2_treat_pileup.bdg -c Peak_Calls/$*/macs2_control_lambda.bdg -o $@ -m logLR -p 0.00001
	sort -k1,1 -k2,2n -o $@ $@

##############
# IDR Metric #
##############
define IDR_Commands

endef

################################
# Peak Call Jaccard Statistics #
################################
define Peak_Heatmap
Results/$1_Peak_Heatmap.png : | $(foreach lib,$($1_LIBRARIES),Peak_Calls/$(lib))
	$(PIPELINE_DIR)/Make_Jaccard_Matrix.sh $1 $2
endef
$(foreach assay,$(NARROW_PEAKS_WITH_INPUT),$(eval $(call Peak_Heatmap,$(assay),narrowPeak)))
$(foreach assay,$(BROAD_PEAKS),$(eval $(call Peak_Heatmap,$(assay),broadPeak)))

define Peak_Summary
Results/$1_Peak_Summary.txt : $(foreach tissue,$($1_TISSUES),Peak_Calls/$1_$(tissue)_Combined_Peaks.bed) | $(foreach lib,$($1_LIBRARIES),Peak_Calls/$(lib))
	python $(PIPELINE_DIR)/Peak_Call_Summary.py $1 "$($1_TISSUES)" "$($1_REPLICATES)" > $$@
endef
$(foreach assay,$(PEAK_ASSAYS),$(eval $(call Peak_Summary,$(assay))))

#######################
# Generate Peak Calls #
#######################
define BroadPeakNoInput
Peak_Calls/$1_% : Aligned_Reads/$1_%.bam | Peak_Calls
	mkdir -p $$@
	macs2 callpeak -t $$< -f BAM -n $$@/macs2 -g $(GENOME_SIZE) -q 0.05 --broad -B --SPMR
	sort -k1,1 -k2,2n $$@/macs2_peaks.broadPeak > $$@/macs2_peaks.bed
endef
$(foreach assay,$(BROAD_PEAKS_NO_INPUT),$(eval $(call BroadPeakNoInput,$(assay))))

define NarrowPeakWithInput
Peak_Calls/$1_% : Aligned_Reads/$1_%.bam Aligned_Reads/Input_%.bam | Peak_Calls
	mkdir -p $$@
	macs2 callpeak -t $$< -c $$(word 2,$$^) -f BAM -n $$@/macs2 -g $(GENOME_SIZE) -q 0.01 -B --SPMR
	sort -k1,1 -k2,2n $$@/macs2_peaks.narrowPeak > $$@/macs2_peaks.bed
endef
$(foreach assay,$(NARROW_PEAKS_WITH_INPUT),$(eval $(call NarrowPeakWithInput,$(assay))))

define BroadPeakWithInput
Peak_Calls/$1_% : Aligned_Reads/$1_%.bam Aligned_Reads/Input_%.bam | Peak_Calls
	mkdir -p $$@
	macs2 callpeak -t $$< -c $$(word 2,$$^) -f BAM -n $$@/macs2 -g $(GENOME_SIZE) -q 0.05 --broad -B --SPMR --fix-bimodal --extsize 200
	sort -k1,1 -k2,2n $$@/macs2_peaks.broadPeak > $$@/macs2_peaks.bed
endef
$(foreach assay,$(BROAD_PEAKS_WITH_INPUT),$(eval $(call BroadPeakWithInput,$(assay))))

######################################
# Long Non-Coding RNA Identification #
######################################
LncRNA/%.fa : LncRNA/%_vcf.gz LncRNA/%_vcf.gz.tbi
	cat LncRNA/All_Transcripts.fa | vcf-consensus $< > $@

LncRNA/%_vcf.gz.tbi : LncRNA/%_vcf.gz
	tabix -p vcf -f $<

LncRNA/%_vcf.gz : LncRNA/%_Sorted.bam
	samtools mpileup -vf LncRNA/All_Transcripts.fa $< > LncRNA/$*_vcf.gz

LncRNA/%_Sorted.bam : LncRNA/%.bam
	samtools sort -@ $(THREADS) -o $@ $<

LncRNA/%.bam : Trimmed_Reads/%_R1_val_1.fq.gz Trimmed_Reads/%_R2_val_2.fq.gz LncRNA/All_Transcripts.1.bt2 
	bowtie2 -x LncRNA/All_Transcripts -1 $< -2 $(word 2,$^) --met-file LncRNA/$*_metrics.txt -p $(THREADS) | samtools view -bS - > $@
	
LncRNA/All_Transcripts.1.bt2 : LncRNA/All_Transcripts.fa
	bowtie2-build $< LncRNA/All_Transcripts

LncRNA/LncRNA_Classes.txt : LncRNA/Intermediate_LncRNA.gtf LncRNA/temp_annotation.gtf
	FEELnc_classifier.pl -i $< -a LncRNA/temp_annotation.gtf > $@

LncRNA/temp_annotation.gtf : $(ANNOTATION)
	python /home/ckern/FAANG_pipelines/gff2gtf.py $< | grep transcript_id | grep protein_coding > $@

LncRNA/Intermediate_LncRNA.gtf : LncRNA/Intermediate_LncRNA_Transcripts.fa Cufflinks_Output/cuffmerge/merged.gtf
	python $(PIPELINE_DIR)/Intersect_Fasta_Gtf.py $^ > $@

LncRNA/Intermediate_LncRNA_Transcripts.fa : LncRNA/Uniprot_Matches.txt LncRNA/Pfam_Matches.txt
	python $(PIPELINE_DIR)/Filter_Fasta.py LncRNA/Uniprot_Matches.txt LncRNA/All_Transcripts.fa > temp && python $(PIPELINE_DIR)/Filter_Fasta.py LncRNA/Pfam_Matches.txt temp > $@

LncRNA/Uniprot_Matches.txt : LncRNA/All_Transcripts.fa
	diamond blastx -p $(THREADS) -q $< -d ~/blastdbs/uniprot_sprot.dmnd -a LncRNA/Uniprot_Matches -e 0.001 -k 1 --sensitive && diamond view -a LncRNA/Uniprot_Matches.daa > $@
	
LncRNA/Pfam_Matches.txt : LncRNA/All_Transcripts.fa
	diamond blastx -p $(THREADS) -q $< -d ~/blastdbs/Pfam-A.dmnd -a LncRNA/Pfam_Matches -e 0.001 -k 1 --sensitive && diamond view -a LncRNA/Pfam_Matches.daa > $@

LncRNA/All_Transcripts.fa : LncRNA/All_Exons.fa
	python $(PIPELINE_DIR)/Merge_Exons.py $< > $@

LncRNA/All_Exons.fa : LncRNA/All_Exons.bed
	bedtools getfasta -fi $(ASSEMBLY) -bed $< -name -s -fo $@

LncRNA/All_Exons.bed : Cufflinks_Output/cuffmerge/merged.gtf | LncRNA
	python $(PIPELINE_DIR)/CuffGTF2Bed.py $< > $@

null :=
space := $(null) #
comma := ,
LncRNA/genes.fpkm : Cufflinks_Output/cuffmerge/merged.gtf | LncRNA
	cuffnorm -o LncRNA -p $(THREADS) -L $(subst $(space),$(comma),$(RNASeq_TISSUES)) -library-type fr-firststrand -library-norm-method geometric $< $(foreach tissue,$(RNASeq_TISSUES),$(subst $(space),$(comma),$(foreach rep,$(RNASeq_$(tissue)_REPLICATES),Aligned_Reads/RNASeq_$(tissue)_$(rep).bam)))


############################
# Cuffmerge Transcriptomes #
############################
Cufflinks_Output/cuffmerge/merged.gtf : Cufflinks_Output/transcriptome_list.txt
	cuffmerge -o Cufflinks_Output/cuffmerge -g $(ANNOTATION) -p $(THREADS) -s $(ASSEMBLY) $<

Cufflinks_Output/transcriptome_list.txt : $(foreach lib,$(RNASeq_LIBRARIES),Cufflinks_Output/$(lib)/transcripts.gtf)
	rm -f $@
	$(foreach lib,$(RNASeq_LIBRARIES),echo Cufflinks_Output/$(lib)/transcripts.gtf >> $@;)


###############################
# Quantify RNA-Seq Expression #
###############################
Bed_Files/%/Expressed_Promoters.bed : Bed_Files/%/Expressed_Genes.bed
	python $(PIPELINE_DIR)/Create_Promoter_Bed.py $< > $@

Bed_Files/%/Repressed_Promoters.bed : Bed_Files/%/Repressed_Genes.bed
	python $(PIPELINE_DIR)/Create_Promoter_Bed.py $< > $@

Bed_Files/%/Expressed_Genes.bed Bed_Files/%/Repressed_Genes.bed : Cufflinks_Output/%/genes.fpkm_tracking | Bed_Files
	mkdir -p Bed_Files/$*
	python $(PIPELINE_DIR)/Make_Gene_Expression_Beds.py $</genes.fpkm_tracking $(ANNOTATION) 1.0 0.2 Bed_Files/$* 	

Cufflinks_Output/%/genes.fpkm_tracking Cufflinks_Output/%/transcripts.gtf : Aligned_Reads/%.bam | Cufflinks_Output
	cufflinks --library-type fr-firststrand -q -L $(subst RNASeq_,,$*) -o Cufflinks_Output/$* -p $(THREADS) -g $(ANNOTATION) $<

############################
# Align Paired-End RNA-Seq #
############################
Aligned_Reads/RNASeq_%.bam : Tophat_Output/RNASeq_%/accepted_hits.bam | Aligned_Reads
	samtools view -b -q 15 $< > $@

Tophat_Output/%/accepted_hits.bam : Trimmed_Reads/%_R1_val_1.fq.gz Trimmed_Reads/%_R2_val_2.fq.gz $(ASSEMBLY_BASE).1.bt2 Temp_Files/transcriptome.1.bt2 | Tophat_Output
	tophat -p $(THREADS) -o Tophat_Output/$* --library-type fr-firststrand --transcriptome-index=Temp_Files/transcriptome $(ASSEMBLY_BASE) $< $(word 2,$^)

Temp_Files/transcriptome.1.bt2 : $(ANNOTATION) $(ASSEMBLY_BASE).1.bt2 | Temp_Files
	tophat -p $(THREADS) -G $(ANNOTATION) --transcriptome-index=Temp_Files/transcriptome $(ASSEMBLY_BASE)

$(ASSEMBLY_BASE).1.bt2 : $(ASSEMBLY)
	bowtie2-build $(ASSEMBLY) $(ASSEMBLY_BASE)

############################################
# Filter ChIP-Seq and DNase-Seq Alignments #
############################################
Aligned_Reads/%.bam : Bwa_Output/%.duplicate-marked.bam | Aligned_Reads
	samtools view -b -q 15 $< > $@

Bwa_Output/%.duplicate-marked.bam : Bwa_Output/%.sorted.bam
	picard-tools MarkDuplicates INPUT=$< OUTPUT=$@ METRICS_FILE=$(@D)/$*.dup_metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT

Bwa_Output/%.sorted.bam : Bwa_Output/%.aligned.bam
	picard-tools SortSam INPUT=$< OUTPUT=$@ SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT

################################
# Align ChIP-Seq and DNase-Seq #
################################
Bwa_Output/%.aligned.bam : Bwa_Output/%_R1_val_1.sai Bwa_Output/%_R2_val_2.sai
	bwa sampe $(ASSEMBLY) $^ Trimmed_Reads/$*_R1_val_1.fq.gz Trimmed_Reads/$*_R2_val_2.fq.gz | samtools view -bS - > $@

Bwa_Output/%.aligned.bam : Bwa_Output/%.sai 
	bwa samse $(ASSEMBLY) $< Trimmed_Reads/$*_trimmed.fq.gz | samtools view -bS - > $@

Bwa_Output/%.sai : Trimmed_Reads/%_trimmed.fq.gz $(ASSEMBLY).bwt | Bwa_Output
	bwa aln -t $(THREADS) -q 15 $(ASSEMBLY) $< > $@

Bwa_Output/%_1.sai : Trimmed_Reads/%_1.fq.gz | Bwa_Output
	bwa aln -t $(THREADS) -q 15 $(ASSEMBLY) $< > $@

Bwa_Output/%_2.sai : Trimmed_Reads/%_2.fq.gz | Bwa_Output
	bwa aln -t $(THREADS) -q 15 $(ASSEMBLY) $< > $@

$(ASSEMBLY).bwt : $(ASSEMBLY)
	bwa index -p $(ASSEMBLY) -a bwtsw $(ASSEMBLY)

##############
# Trim Reads #
##############
Trimmed_Reads/%_R1_val_1.fq.gz Trimmed_Reads/%_R2_val_2.fq.gz : Raw_Reads/%_R1.fq.gz Raw_Reads/%_R2.fq.gz | Trimmed_Reads
	trim_galore --paired $< $(word 2,$^) -o Trimmed_Reads

Trimmed_Reads/%_trimmed.fq.gz : Raw_Reads/%.fq.gz | Trimmed_Reads
	trim_galore $< -o Trimmed_Reads

# Create missing directories
$(DIRS) : 
	mkdir $@

%.bam.bai : %.bam
	samtools index $<

clean :
	rm -rf Bwa_Output Trimmed_Reads Tophat_Output DeepTools Bed_Files

deepclean : clean
	rm -rf Trimmed_Reads Aligned_Reads Peak_Calls

.PHONY : $(PHONIES) clean deepclean
