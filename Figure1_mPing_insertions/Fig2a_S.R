pdf("Fig2a_S.pdf")
x <- read.table("RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.type.summary2")
pie(x$V2, labels=x$V1, radius=1, col=c("purple", "blue", "red", "green"))
dev.off()
