pdf("Excision_distance.matrix_events.plot.pdf")
par(mar=c(5,5,5,2))
read.table("Excision_distance_RIL.matrix_events.1.txt", header=T) -> x
x <- x[x[,3]<1200000,]
plot(x[,3]/1000, x[,2], type="p", pch=18, col="black", xlab="Distance (kb)", ylab="Number of excision", xlim=c(min(x[,3])/1000, 1200000/1000), ylim=c(0, 20), cex.axis=1.4, cex.lab=1.4, axes=FALSE)
axis(2, at=seq(0, 20, by=5),labels=seq(0, 20, by=5), cex=1.4, cex.axis=1.4)
axis(1, at=c(round(min(x[,3])/1000, 2), seq(200, 1200, by=200)),labels=c(round(min(x[,3])/1000, 2), seq(200, 1200, by=200)), cex=1.4, cex.axis=1.4)
x <- x[x[,3]<100000,]
plot(x[,3]/1000, x[,2], type="p", pch=18, col="black", xlab="Distance (kb)", ylab="Number of excision", xlim=c(min(x[,3])/1000, 100000/1000), ylim=c(0, 20), cex.axis=2, cex.lab=2, axes=FALSE)
axis(2, at=seq(0, 20, by=5),labels=seq(0, 20, by=5), cex=2, cex.axis=2, cex.lab=2)
axis(1, at=c(round(min(x[,3])/1000, 2), seq(20, 100, by=20)),labels=c(round(min(x[,3])/1000, 2), seq(20, 100, by=20)), cex=2, cex.axis=2)
dev.off()

