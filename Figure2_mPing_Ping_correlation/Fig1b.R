error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y, angle=90, code=1, length=length, ...)
    #text(x, y+upper, paste("n=", ))
}

pdf("Fig1b.pdf")
par(mar=c(5,5,4,2))
data =read.table("RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.ping_number.summary.1")
expr = data[,2]
std = data[,3]
barx <- barplot(expr, col=c("cornflowerblue"), ylim=c(0,120), border=F, axis.lty=1, xlab='', ylab='')
error.bar(barx, expr, std)
axis(1,c(0.1, max(barx)+0.6),line=0,labels=c("",""), cex=1.4)
text(barx, rep(-5, 6),offset=2,labels=data[,1],srt=0,xpd=TRUE, cex=1.4)
#legend("topright",c("HEG4","Nipponbare"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("blue","orange"))
xpos <- 3.6
ypos <- 44
mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1.4, col="black")
mtext("copy number", side=1,font=1, at=xpos+2.1,line=3, cex=1.4, col="black")
mtext("Unique homozygous", side=2,font=1, at=ypos,line=3, cex=1.4, col="black")
mtext("mPing", side=2,font=3, at=ypos+34,line=3, cex=1.4, col="black")
mtext("number", side=2,font=1, at=ypos+53,line=3, cex=1.4, col="black")
#text(2, 115, 'Pearson correlation:', cex=1.4)
#text(3.7, 108, 'R-squared = 0.70, p-value < 2.2e-16', cex=1.4)
#n
text(barx, expr+std+10, paste('n = ', data[,6], sep=''), cex=1.4)
dev.off()


