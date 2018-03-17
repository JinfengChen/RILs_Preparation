pdf('mping.excision_events.distr.pdf')
t <- read.table("Excision_newpipe_version1.footprint.list.draw.txt")
f <- t[,2]
f <- f[f>0]
brk <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
xx0 <- hist(f,breaks=brk, ylab="Frequency", xlab="Excisions", col="steelblue2", main="",border=T, plot=F)
y <- xx0$counts
x <- xx0$mids

xx <- barplot(y,beside=TRUE,ylab="Number of mPing loci (n=178)",cex.lab=1.2, cex.names=1.2, cex.axis=1.2,border=FALSE,axes=FALSE,ylim=c(0,120),col="steelblue2")

b <- 220

text(seq(0.7,22, by=1.2), rep(-5, 18), labels=c(1:18), xpd=TRUE, cex=1.2)
axis(1, c(0,21.7), line=0, labels=c("",""))
axis(2, c(0,20,40,60,80,100,120), line=0, labels=c(0,20,40,60,80,100,120))
rect(0.2,b,max(xx)+0.8,b+5,border=FALSE,col='white')
#text(c(0,2.5,4.8,7.3,9.7,12.1,14.5,16.9),rep(-45,length(y)), cex=1, offset=2,labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),srt=0,xpd=TRUE)
mtext("Number of excision events", side=1,cex=1.2, at=10,line=3)
#n
text(xx, y+3, y, cex=1.2, xpd=TRUE)
#cutoff
abline(v=4.9, col='gray', lty=2, lwd=2)
#p-value
text(5.6, 116, 'Binomial test, P-value=6.31e-5', adj=0, cex=1.2)
dev.off()
