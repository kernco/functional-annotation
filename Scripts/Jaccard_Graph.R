library(gplots)
library(RColorBrewer)

jaccard_table = read.table(snakemake@input[[1]])
jaccard_matrix = as.matrix(jaccard_table)
png(snakemake@output[[1]])
heatmap.2(jaccard_matrix, col=colorRampPalette(brewer.pal(9,"Reds")), margins=c(14,14), density.info="none", main=strsplit(snakemake@wildcards[["assay"]],'.',fixed=TRUE)[[1]][1], lhei=c(2,8), trace="none", dendrogram='column')
dev.off()
