library(normr)
#library("IRanges")
genome <- read.table(snakemake@input[["lengths"]])
#gr <- GenomicRanges::GRanges(as.character(genome[,1]), IRanges(1,as.integer(genome[,2])))
seqlengths <- genome[,2]
names(seqlengths) <- genome[,1]
gr <- GenomicRanges::tileGenome(seqlengths, tilewidth=200, cut.last.tile.in.chrom=T)
print(length(gr))
chipcounts <- read.table(snakemake@input[["treatment"]], sep="\t")[[4]]
print(length(chipcounts))
contcounts <- read.table(snakemake@input[["control"]], sep="\t")[[4]]
fit <- regimeR(treatment=chipcounts, control=contcounts, genome=gr, models=3, verbose=FALSE)
summary(fit)
exportR(fit, filename=snakemake@output[["regions"]], type="bed", fdr=snakemake@params[["fdr"]])
