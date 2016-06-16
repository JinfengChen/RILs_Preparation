pdf("NB_Ping_scatter_plot.pdf")
par(mar=c(5,5,4,2))
mping <- read.table("None.table.txt", header=TRUE)
beeswarm(Heterozygous_mPing ~ Code, data = mping,
          pch = 16, col=c("cornflowerblue"),
          xlab = "", ylab = "",
          labels = c(2,3,4,5,6,7), axes=FALSE)
bxplot(Heterozygous_mPing ~ Code, data = mping, add=TRUE)
#axis(2,seq(0, 140, by=20),line=0, labels=seq(0, 140, by=20), cex.axis=1.4)
#axis(1,c(0.5, 8.5),line=0,labels=c("",""), cex.axis=1.4)
#x and y lab
xpos=3
ypos=50
mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1.4, col="black")
mtext("copy number", side=1,font=1, at=xpos+1.8,line=3, cex=1.4, col="black")
mtext("Unique homozygous", side=2,font=1, at=ypos,line=3, cex=1.4, col="black")
mtext("mPing", side=2,font=3, at=ypos+44,line=3, cex=1.4, col="black")
mtext("number", side=2,font=1, at=ypos+68,line=3, cex=1.4, col="black")
#x annotation
text(c(1,2,3,4,5,6,7,8), rep(-15, 7), offset=2,labels=c(0,1,2,3,4,5,6,7),srt=0,xpd=TRUE, cex=1.4)

dev.off()
