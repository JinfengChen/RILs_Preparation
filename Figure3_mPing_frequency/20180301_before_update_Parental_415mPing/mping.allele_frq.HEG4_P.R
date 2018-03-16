pdf('mping.allele_frq.HEG4_P.pdf')
t <- read.table("HEG4.2.3.RelocaTE2.mping.non-ref.allele.frq.sorted")
f <- t[,3]
#f <- f[f>0]
breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8)
xx <- hist(f, breaks=breaks,ylab="mPing insertions", xlab="Allele Frequency", col="steelblue2", main="",border=T, plot=F)
y <- xx$counts
x <- xx$mids

xx <- barplot(y,beside=TRUE, ylab="Number of mPing loci",cex.lab=1.2, cex.names=1.2, cex.axis=1.2,border=FALSE,axes=FALSE, ylim=c(0, 250), col="steelblue2")

b <- 220

axis(1,c(0,2.5,4.8,7.3,9.7,12,14.5,17, 19.3),line=0,c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))
axis(2,c(0,50,100,150,200,250),line=0,labels=c(0,50,100,150,200,250))
#rect(0.2,b,max(xx)+0.8,b+5,border=FALSE,col='white')
#text(c(0,2.5,4.8,7.3,9.7,12.1,14.5,16.9),rep(-45,length(y)), cex=1, offset=2,labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),srt=0,xpd=TRUE)
mtext("HEG4 mPing allele frequency", side=1,cex=1.2, at=9,line=3)

dev.off()
