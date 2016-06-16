pdf("Fig1b_unique_mping_hist.pdf")

par(mar=c(5,5,4,2))
mping <- read.table("RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt", header=TRUE)
hist(mping$Unique_hom, breaks=20, xlab="", ylab="Number of RILs", col="cornflowerblue", main="", cex.axis=1.4, cex.lab=1.4)
xpos=40
mtext("Unique homozygous", side=1,font=1, at=xpos,line=3, cex=1.4, col="black")
mtext("mPing", side=1,font=3, at=xpos+40,line=3, cex=1.4, col="black")
mtext("number", side=1,font=1, at=xpos+63,line=3, cex=1.4, col="black")

dev.off()
