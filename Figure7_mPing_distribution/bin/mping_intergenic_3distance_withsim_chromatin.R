error.bar <- function(x, y, upper, lower=upper, color,length=0.06,...){
     if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
     stop("vectors must be same length")
     arrows(x,y+upper, x, y-lower, col=color,angle=90, code=3, length=length, ...)
 }

pdf("mping_intergenic_3distance_withsim_chromatin.pdf")

par(mar=c(6,4,4,2), cex=1.2)
ril5 <- read.table("../../Figure7_mPing_distribution_pre/bin/RIL.mRNA.3primer.distance.distr")
sim_both5 <- read.table("../../Figure7_mPing_distribution_pre/bin/results_simulationV2_Chromatin0.99_TSD9mer_ril/Simulate.TSD9mer.rilMat.mRNA.3primer.distance.distr")
sim_dhs5 <- read.table("../../Figure7_mPing_distribution_pre/bin/results_simulationV2_Chromatin0.99_ril/Simulate.TSD9mer.rilMat.mRNA.3primer.distance.distr")
sim_tsd5 <- read.table("../../Figure7_mPing_distribution_pre/bin/results_simulationV2_TSD9mer_ril/Simulate.TSD9mer.rilMat.mRNA.3primer.distance.distr")
control5 <- read.table("../../Figure7_mPing_distribution_pre/bin/results_simulationV2_Random_TSD9mer_ril/Simulate.TSD9mer.rilMat.mRNA.3primer.distance.distr")

ril5 <- subset(ril5, V1 >= -3 & V1 <= 4)
sim_both5 <- subset(sim_both5, V1 >= -3 & V1 <= 4)
sim_dhs5  <- subset(sim_dhs5, V1 >= -3 & V1 <= 4)
sim_tsd5  <- subset(sim_tsd5, V1 >= -3 & V1 <= 4)
control5  <- subset(control5, V1 >= -3 & V1 <= 4)

#ril5 <- ril5[-length(ril5[,1]),]
#sim_both5 <- sim_both5[-length(sim_both5[,1]),]
#sim_dhs5 <- sim_dhs5[-length(sim_dhs5[,1]),]
#sim_tsd5 <- sim_tsd5[-length(sim_tsd5[,1]),]
#control5 <- control5[-length(control5[,1]),]

x <- c(-1750, -1250, -750, -250, 250, 750, 1250, 1750)
plot(x, ril5[,4], type='b', pch= 1,lwd = 2 , col="aquamarine3", xaxt='n', frame.plot = FALSE, xlim=c(-2000, 2000), ylim=c(0,0.2), ylab="Proportion", xlab="")
lines(x, sim_both5[,4], type='b',pch= 2,lwd = 2 , col="steelblue2")
error.bar(x, sim_both5[,4], sim_both5[,7]-sim_both5[,4], sim_both5[,7]-sim_both5[,4], 'steelblue2')

lines(x, sim_dhs5[,4], type='b',pch= 3,lwd = 2 , col="sandybrown")
error.bar(x, sim_dhs5[,4], sim_dhs5[,7]-sim_dhs5[,4], sim_dhs5[,7]-sim_dhs5[,4], 'sandybrown')

lines(x, sim_tsd5[,4], type='b',pch= 20, cex=0.2,lwd = 2 , col="darkred")
error.bar(x, sim_tsd5[,4], sim_tsd5[,7]-sim_tsd5[,4], sim_tsd5[,7]-sim_tsd5[,4], 'darkred')

lines(x, control5[,4], type='b',pch= 20, cex=0.2,lwd = 2 , col="dim gray")
error.bar(x, control5[,4], control5[,7]-control5[,4], control5[,7]-control5[,4], 'dim gray')

#axis(1,seq(1:length(ril5[,1])),line=0, labels=rep("",length(ril5[,1])))
#text(seq(1:length(ril5[,1][-1]))+0.5,rep(-0.02,7), cex=1, offset=2,labels=ril5[,1][-length(ril5[,1])]*500/-1000,srt=55,xpd=TRUE)

axis(1,seq(-2000, 2000, by= 500),line=0, labels=seq(-2, 2, by= 0.5))
legend('topright', bty='n', border='NA', lty= c(1,2,3,4,5), pch = c(1,2,3,4,20), cex=1 , lwd = 2 ,col=c("aquamarine3", "steelblue2", "sandybrown", "darkred", "dim gray"), c("Unique", "Sim_DHS_TSD", "Sim_DHS", "Sim_TSD", "Control"))
mtext("Distance to TTS (kp)", side=1,cex=1.2, at=9,line=3)



dev.off()

