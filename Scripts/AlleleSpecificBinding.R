library(BaalChIP)

samplesheet <- snakemake@input[["samplesheet"]]
hets <- c("GRP"=snakemake@input[["hets"]])
res <- BaalChIP(samplesheet=samplesheet, hets=hets, CorrectWithgDNA=snakemake@input[["gdna"]])
res <- alleleCounts(res)
#Possibly add filtering steps here
res <- mergePerGroup(res)
res <- filter1allele(res)
res <- getASB(res, Iter=5000, conf_level=0.95, cores=4, RMcorrection=FALSE, RAFcorrection=FALSE)
result <- BaalChIP.report(res)
qc <- summaryQC(res)[["filtering_stats"]]
asb <- summaryASB(res)
write.table(result[result$isASB==TRUE,], file=snakemake@output[["table"]])
write.table(qc, file=snakemake@output[["qc"]])
write.table(asb, file=snakemake@output[["summary"]])
