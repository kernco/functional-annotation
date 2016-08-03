#!/bin/bash -l
#
# Usage: sbatch [sbatch-options] slurm_script.sh Makefile_target Jobs Threads
# Jobs is the number of jobs make should run in parallel
# Threads is the number of threads each job can use
# Make sure you have Jobs x Threads CPU cores available
# Too many jobs in parallel may hinder performance due to excessive disk access

#SBATCH -J FAANG-Pipeline

module load trim-galore tophat cufflinks bowtie2 samtools bwa python pysam picardtools deepTools R numpy scipy matplotlib vcftools diamond
PYTHONPATH=$PYTHONPATH:/home/ckern/deepTools
FEELNCPATH=/home/ckern/FEELnc
PERL5LIB=$PERL5LIB:${FEELNCPATH}/lib/
PATH=$PATH:${FEELNCPATH}/scripts/:${FEELNCPATH}/bin/LINUX/

make -k -j $2 -f ~/Integrative_Pipeline/Makefile $1 THREADS=$3
