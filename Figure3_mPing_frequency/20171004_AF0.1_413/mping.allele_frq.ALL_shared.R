library("plotrix")
pdf('mping.allele_frq.ALL_shared.pdf')
t <- read.table("RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency")
f <- t[,7]
#f <- f[f>0.005]
breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8)
xx <- hist(f,breaks=breaks,ylab="Insertion sites",xlab="Allele Frequency", col="steelblue2", main="",border=T, plot=F)
y <- xx$counts
x <- xx$mids
cut <- 1500
y[1] <- y[1] - cut

xx <- barplot(y,beside=TRUE,ylab="Number of mPing loci", cex.lab=1.2, cex.names=1.2, cex.axis=1.2, border=FALSE,axes=FALSE, ylim=c(0,400), col="steelblue2")

b <- 240

axis(1,c(0,2.5,4.8,7.3,9.7,12,14.5,17),line=0,c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
axis(2,c(0,100,200,300,400),line=0,labels=c(0,100,200,1800,1900))
axis.break(2, b,style="slash")
rect(0.2,b,max(xx)+0.8,b+5,border=FALSE,col='white')
#text(c(0,2.5,4.8,7.3,9.7,12.1,14.5,16.9),rep(-45,length(y)), cex=1, offset=2,labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),srt=0,xpd=TRUE)
mtext("Shared mPing allele frequency in RILs", side=1,cex=1.2, at=9,line=3)

dev.off()
