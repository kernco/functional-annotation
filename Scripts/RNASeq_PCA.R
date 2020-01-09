library(ggfortify)

rep1 <- read.table(snakemake@input['reps'][0], sep="\t")
rep2 <- read.table(snakemake@input['reps'][1], sep="\t")
table <- cbind(rep1, rep2)
pc <- prcomp(test, scale = T, center=T)
autoplot(pc)
ggsave(snakemake@output[png])
