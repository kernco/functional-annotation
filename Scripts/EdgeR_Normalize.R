library(edgeR)

RG <- readDGE(snakemake@input[["htseqs"]], group=snakemake@params[["groups"]], labels=snakemake@params[["labels"]])
metatags <- grep("^__", rownames(RG))
RG <- RG[-metatags, , keep.lib.sizes=FALSE]
keep <- filterByExpr(RG)
RG <- RG[keep, , keep.lib.sizes=FALSE]
RG <- calcNormFactors(RG)
RG <- cpm(RG)
for(rep in 1:length(snakemake@params[["replicates"]])) {
    write.table(RG[,grepl(paste0("_",snakemake@params[["replicates"]][rep],"$"), colnames(RG))], file=snakemake@output[[rep]], sep="\t", quote=F)
}
