pdf("TSSprofile_Nuc_DHS.pdf", height=6, width=8)
par(mar=c(6,4,4,2), cex=1.2)
nuc <- read.table("Nucleosome_TSS.profile")
dhs <- read.table("DHS_TSS.profile")
#ngene <- 5157  #chr1
#nreads_nuc <- 7480914/1000000
#nreads_dhs <- 3279088/1000000
bandwidth_nuc  <- 20 #bandwidth used to smooth the curve
bandwidth_dhs  <- 150

ngene <- 39900 #all
nreads_nuc <- 60745783/1000000
nreads_dhs <- 23299296/1000000
#normalized read count per base per million reads
#Nucleosome
plot(nuc[,1]/nreads_nuc/ngene/bandwidth_nuc, type='l', pch= 1,lwd = 2 , col="brown", xaxt='n', frame.plot = FALSE, ylim=c(0.001, 0.004), ylab="Normalized MNase-seq Reads", xlab="")
axis(1,seq(0, 4000, by= 500),line=0, labels=seq(-2, 2, by= 0.5))
mtext("Distance to TSS (kb)", side=1, cex=1.2, at=2000, line=3)
#DHS
plot(dhs[,1]/nreads_dhs/ngene/bandwidth_dhs, type='l', pch= 1,lwd = 2 , col="blue", xaxt='n', frame.plot = FALSE, ylim=c(0.002, 0.006), ylab="Normalized DNase-seq Reads", xlab="")
axis(1,seq(0, 4000, by= 500),line=0, labels=seq(-2, 2, by= 0.5))
mtext("Distance to TSS (kb)", side=1, cex=1.2, at=2000, line=3)
#combined
plot(nuc[,1]/nreads_nuc/ngene/bandwidth_nuc, type='l', pch= 1,lwd = 2 , col="brown", xaxt='n', frame.plot = FALSE, ylim=c(0.001, 0.006), ylab="Normalized Reads count", xlab="")
lines(dhs[,1]/nreads_dhs/ngene/bandwidth_dhs, type='l', pch= 1, lwd = 2 , col="blue")
axis(1,seq(0, 4000, by= 500),line=0, labels=seq(-2, 2, by= 0.5))
#axis(4,seq(0.01, 0.))
#text(seq(1:length(dist[,1])),rep(-0.04,7), cex=1, offset=2,labels=rev(dist[,1]/1000),srt=0,xpd=TRUE)
mtext("Distance to TSS (kb)", side=1, cex=1.2, at=2000, line=3)
legend('topright', bty='n', border='NA', lty= c(1,1), cex=1 , lwd = 2 ,col=c("brown", "blue"), c("MNase-seq reads", "DNase-seq reads"))

dev.off()

