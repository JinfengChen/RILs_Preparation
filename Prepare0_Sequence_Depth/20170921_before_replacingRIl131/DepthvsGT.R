pdf("DepthvsGT.pdf", height=7, width=7)
par(mar=c(5,5,4,2))
read.table("Depth.table", header=TRUE) -> x
x <- x[x[,3] < 30,]
plot(x[,3],1-x[,4],type="p",pch=18,col="black",xlab="Read Depth (X)",ylab="Genotyped SNPs (%)",xlim=c(0,30),ylim=c(0,1), axes=FALSE, cex.lab=1.4, cex=1.4)
atx <- seq(0,30, by=5)
aty <- seq(0,1, by=0.2)
axis(1,at=atx,labels=atx, cex.axis=1.4, cex=1.4)
axis(2,at=aty,labels=aty, cex.axis=1.4, cex=1.4)
dev.off()

