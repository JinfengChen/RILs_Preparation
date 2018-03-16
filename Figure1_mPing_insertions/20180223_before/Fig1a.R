pdf("Fig1a.pdf", height=7, width=6)
par(mfrow=c(1,2))

##type of mping insertion, parental, shared, unique
x <- read.table("RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.type.summary1", row.names=1)
x <- t(x)
data <- rbind(x[,1]/sum(x[1,]), x[,2]/sum(x[1,]))
rownames(data) <- colnames(x)

#labels percentage and number: paste(x[1], "(", round(data[1,1]*100), "%)", sep="")
par(mar=c(5,5,4,5))
barplot(data, ylab="Proportion", space=0.1, col=c("darkseagreen", "lightpink3"), cex.axis=1.2, cex.names=1.2, cex.lab=1.2)
legend(-0.1, 1.17,c("Shared", "Unique"), bty="n", xpd=TRUE, border="NA",lty=c(0,0), cex=1.2, fill=c("lightpink3", "darkseagreen"))
text(0.6, 0.355, cex=1.4, paste("(", x[1], ")", sep=""), col='white')
text(0.6, 0.4, cex=1.4, paste(round(data[1,1]*100), "%", sep=""), xpd=TRUE, col='white')
text(0.6, 0.895, cex=1.4, paste("(", x[2], ")", sep=""), col='white')
text(0.6, 0.94, cex=1.4, paste(round(data[2,1]*100), "%", sep=""), xpd=TRUE, col='white')
#text(1.5, 0.96, cex=1.4, paste("(", x[3], ")", sep=""), col='black', xpd=TRUE)
#text(1.5, 1, cex=1.4, paste(round(data[3,1]*100), "%", sep=""), xpd=TRUE, col='black')
text(0.3, -0.08, cex=1.2, 'de novo mPing', font=3, xpd=TRUE)
text(1.51, -0.075, cex=1.2, 'insertions', font=1, xpd=TRUE)
text(0.6, -0.13, cex=1.2, '(n = 16,452)', xpd=TRUE)

##class of mping insertion, homozygous, heterozygous
x <- read.table("RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.class.summary", row.names=1, header=T)
x[,2] <- x[,2] + x[,3]
x[,3] <- NULL
x <- t(x)
data <- cbind(x[,1]/sum(x[,1]), x[,2]/sum(x[,2]), x[,3]/sum(x[,3]))
colnames(data) <- colnames(x)

data1 <- data[,3]
y <- rbind(data1[1], data1[2])
rownames(y) <- c("Homozygous", "Heterzygous")
colnames(y) <- ""
#percentage and number
barplot(y, xlab="", ylab="", col=c("cornflowerblue", "bisque4"), axes=FALSE, cex.axis=1.4, cex.lab=1.4)
legend(0, 1.15,c("Homozygous", "Heterozygous"),bty="n", xpd=TRUE, border="NA",lty=c(0,0),cex=1.2,fill=c("cornflowerblue", "bisque4"))
#text(0.7, 0.4, cex=0.8, paste(x[1,3], "(", round(y[1]*100), "%)", sep=""))
#text(0.7, 0.9, cex=0.8, paste(x[2,3], "(", round(y[2]*100), "%)", sep=""))
#text(0.7, 0.4, cex=1.4, x[1,3], col='white')
#text(1.5, 0.4, cex=1.4, paste(round(y[1]*100), "%", sep=""), xpd=TRUE)
#text(0.7, 0.9, cex=1.4, x[2,3], col='white')
#text(1.5, 0.9, cex=1.4, paste(round(y[2]*100), "%", sep=""), xpd=TRUE)
text(0.7, 0.355, cex=1.4, paste("(", x[1,3], ")", sep=""), col='white')
text(0.7, 0.4, cex=1.4, paste(round(y[1]*100), "%", sep=""), xpd=TRUE, col='white')
text(0.7, 0.855, cex=1.4, paste("(", x[2,3], ")", sep=""), col='white')
text(0.7, 0.9, cex=1.4, paste(round(y[2]*100), "%", sep=""), xpd=TRUE, col='white')
text(0.7, -0.075, cex=1.2, 'Unique', xpd=TRUE)

dev.off()

