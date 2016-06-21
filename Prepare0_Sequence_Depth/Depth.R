pdf("Depth.pdf", height=7, width=7)
par(mar=c(5,5,4,2))
read.table("Depth.table", header =TRUE) -> x
br = c(seq(0, 28, by = 3), 100)
hist(x[,2], breaks=br, plot = FALSE) -> xh
hist(x[,3], breaks=br, plot = FALSE) -> xh1
xx <- xh$mids
xx[length(xx)] <- 28.5
atxvalue <- seq (0, 31,by=3)
atyvalue <- seq (0, 100, by=20)
atyvalue1 <- seq (0, 100, by=20)
atx <- atxvalue
atx_label <- atx
#atx_label[length(atx_label)] <- '>30'
aty <- c(atyvalue, 100)
plot(xx,xh$counts,type="b",pch=18,col="lightseagreen",xlab="Read Depth (X)",ylab="Number of RILs",xlim=c(0,30),ylim=c(0, 100),axes=FALSE, cex.lab=1.4)
lines(xx,xh1$counts,type="b",pch=19,col="tomato4")
legend(14,80,c("Raw Depth", "Mapped Depth"), pch=c(18,19),lty=c(1,1), col=c("lightseagreen", "tomato4"), border=FALSE,bty = 'n', cex=1.4)
axis(1,at=atx,labels=atx_label, cex=1.4, cex.axis=1.4)
axis(2,at=aty,labels=aty, cex=1.4, cex.axis=1.4)
dev.off()

