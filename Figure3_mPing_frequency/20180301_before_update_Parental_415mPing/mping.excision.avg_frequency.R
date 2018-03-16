error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y,angle=90, code=1, length=length, ...)
}

library("plotrix")
pdf('mping.excision.avg_frequency.pdf')
t <- read.table("mping.excision.avg_frequency.table")
std <- t[,3][1:14]
y <- t[,2][1:14]
#y <- append(y,0,0)
#std <- append(std,0,0)
x <- t[,1]

xx <- barplot(y,beside=TRUE,ylab="Number of excised RILs", cex.lab=1.2, cex.names=1.2,cex.axis=1.2,border=FALSE,axes=FALSE, ylim=c(0,80),col="steelblue2")

error.bar(xx, y, std)

b <- 220

axis(1,c(0,2.5,4.8,7.3,9.7,12,14.5,17),line=0,c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))
#axis(1,c(0,2.5,4.8,7.3,9.7,12,14.5,17),line=0,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7))
axis(2,c(0,20,40,60, 80),line=0,labels=c(0,20,40,60, 80))
rect(0.2,b,max(xx)+0.8,b+5,border=FALSE,col='white')
#text(c(0,2.5,4.8,7.3,9.7,12.1,14.5,16.9),rep(-45,length(y)), cex=1, offset=2,labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),srt=0,xpd=TRUE)
mtext("Parental mPing allele frequency", side=1,cex=1.2, at=9,line=3)

dev.off()
