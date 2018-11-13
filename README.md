# UC Davis FAANG Functional Annotation Pipeline
Pipeline for functional annotation of genomes by integrating data from RNA-seq, ChIP-seq, DNase-seq, etc.

This is a work-in-progress. Documentation is currently very brief and minimal. 

## Setting up your data

Create a directory that will be the working directory for the pipeline. You will move your
raw data into this folder, and the pipeline will populate it with many other folders containing
results, analyses, etc.

Inside your working directory, create a folder named `Raw_Reads` and move your raw FastQ files
into this folder. The pipeline gets a lot of information from the filename of the FastQ files,
so they should be named in a specific manner. Each filename should follow the pattern
`Assay_Group_Replicate.fq.gz` where Assay denotes the type of data such as CTCF, H3K4me3, or RNASeq.
The Assay can be named whatever you wish (you will define how each assay should be analyzed
later in a configuration file), with the exception of RNA-seq data which in the current version
of the pipeline must be named RNASeq. The Group part of the filename should be a unique identifier
that all biological replicates share. This is usually something like tissue, breed, or condition.
Finally the Replicate portion of the filename is an identifier for the organism the sample from.
Here are some examples of raw data filenames:

```
H3K4me3_Liver_A.fq.gz
CTCF_LeghornTreated_1163.fq.gz
```

## Configuration file

In your working directory should be a file named config.yaml. An example config.yaml is included
with this pipeline on Github, which can be copied into your working directory and modified
for your project. More detailed documentation about the config.yaml file will be added in the
future, but hopefully the example file will be self-explanatory enough to get it up and
running with your data.

## Running the pipeline

The pipeline can be run with this command in your working directory:

`snakemake -j 4 -s /path/to/pipeline/Snakefile --configfile config.yaml --use-conda`

The -j arugment specifies how many tasks to run in parallel. The -s argument should be
the path to where the Snakefile file of the pipeline is on your system.

If you are running the pipeline on a computing cluster, you can use the run.sh script
included with the pipeline:

`/path/to/pipeline/run.sh`

This is configured for a cluster using the SLURM job scheduler. If your cluster uses
another job scheduler you will need to modify this script. The script itself can be
submitted as a job, which will then submit additional jobs for each step of the pipeline.

## Results

The folders `Figures` and `Tables` will contain many results of the pipeline. BAM files for
the read alignments are in `Aligned_Reads` and peak call files are contained in `Peak_Calls`.

Currently the ChromHMM portion of the pipeline does not use a conda environment and includes
paths to where ChromHMM is installed on our cluster at UC Davis. It is possible to modify
the pipeline for where ChromHMM is on your own system and get it to work, but in the future
this will be made much easier.


