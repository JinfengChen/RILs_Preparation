
pdf("Ping4_hom.pdf")
par(mar=c(5,5,4,2))
library(beeswarm)
mping <- read.table("Ping4_hom.table.txt", header=TRUE)
beeswarm(Heterozygous_mPing ~ Code, data = mping,
          pch = 16, col=c("cornflowerblue"),
          xlab = "", ylab = "",
          labels = c(2,3,4,5,6,7), axes=FALSE, ylim=c(0, 100), xlim=c(0.5, 12.5))
bxplot(Heterozygous_mPing ~ Code, data = mping, add=TRUE)
axis(2,seq(0, 100, by=20),line=0, labels=seq(0, 100, by=20), cex.axis=1.4)
axis(1,c(0.5, 12.5),line=0,labels=c("",""), cex.axis=1.4)
#x and y lab
xpos=5
ypos=30
#mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1.4, col="black")
#mtext("copy number", side=1,font=1, at=xpos+2.6,line=3, cex=1.4, col="black")
mtext("Unique homozygous", side=2,font=1, at=ypos,line=3, cex=1.4, col="black")
mtext("mPing", side=2,font=3, at=ypos+31,line=3, cex=1.4, col="black")
mtext("number", side=2,font=1, at=ypos+48,line=3, cex=1.4, col="black")
#x annotation
text(c(1,2,3,4,5,6,7,8,9,10,11,12), rep(-11, 12), offset=2,labels=c("+", "-", "+", "-","+", "-","+", "-","+", "-", "+", "-","+", "-","+", "-"),srt=0,xpd=TRUE, cex=1.4)
text(c(1.5,3.5,5.5,7.5,9.5,11.5), rep(-19, 12), offset=2,labels=c("2 Pings", "3 Pings", "4 Pings", "5 Pings", "6 Pings", "7 Pings"),srt=0,xpd=TRUE, cex=1.4)

dev.off()

