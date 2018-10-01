#!/bin/bash -l

module load bio3

mkdir -p Logs
snakemake -j 24 --cluster-config /home/ckern/functional-annotation/cluster.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} -n {cluster.cpus} -J {cluster.name} -o {cluster.output} -e {cluster.output}" -s /home/ckern/functional-annotation/Snakefile --configfile config.yaml --latency-wait 120 -p -k --use-conda $@

