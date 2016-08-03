library(gplots)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1])
jaccard_table = x[, -1]
jaccard_matrix = as.matrix(jaccard_table)
png(args[2])
heatmap.2(jaccard_matrix, col=colorRampPalette(brewer.pal(9,"Reds")), margins=c(14,14), density.info="none", main=strsplit(args[1],'.',fixed=TRUE)[[1]][1], lhei=c(2,8), trace="none")
dev.off()
