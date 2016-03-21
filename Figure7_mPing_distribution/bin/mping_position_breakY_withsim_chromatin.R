error.bar <- function(x, y, upper, lower=upper, color,length=0.06,...){
     if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
     stop("vectors must be same length")
     arrows(x,y+upper, x, y-lower, col=color,angle=90, code=3, length=length, ...)
 }


pdf("mping_position_breakY_withsim_chromatin.pdf")
library("plotrix")
par(mar=c(6,4,4,2), cex=1.2)

ril5 <- read.table("../../Figure7_mPing_distribution_pre/bin/RIL.position.distr")
sim_both5 <- read.table("../../Figure7_mPing_distribution_pre/bin/results_simulationV2_Chromatin0.99_TSD9mer_ril/Simulate.TSD9mer.rilMat.position.distr")
sim_dhs5 <- read.table("../../Figure7_mPing_distribution_pre/bin/results_simulationV2_Chromatin0.99_ril/Simulate.TSD9mer.rilMat.position.distr")
sim_tsd5 <- read.table("../../Figure7_mPing_distribution_pre/bin/results_simulationV2_TSD9mer_ril/Simulate.TSD9mer.rilMat.position.distr")
control5 <- read.table("../../Figure7_mPing_distribution_pre/bin/results_simulationV2_Random_TSD9mer_ril/Simulate.TSD9mer.rilMat.position.distr")

ril5p <- ril5[,3]/sum(ril5[,3])
sim_both5p <- sim_both5[,3]/sum(sim_both5[,3])
sim_dhs5p <- sim_dhs5[,3]/sum(sim_dhs5[,3])
sim_tsd5p <- sim_tsd5[,3]/sum(sim_tsd5[,3])
control5p <- control5[,3]/sum(control5[,3])

#cut <- 0.5
#x1p[1] <- x1p[1] - cut
#x2p[1] <- x2p[1] - cut
#x3p[1] <- x3p[1] - cut
#x4p[1] <- x4p[1] - cut

data <- rbind(ril5p, sim_both5p, sim_dhs5p, sim_tsd5p, control5p)
xx <- barplot(data,beside=TRUE,ylab="Proportion",cex.names=1,cex.axis=1,border=FALSE,axes=FALSE,ylim=c(0,0.5),col=c("aquamarine3", "steelblue2" ,"sandybrown", "darkred", "dim gray"))
#error.bar(xx[4,], x4p, x4[,7]-x4[,4], x4[,7]-x4[,4], 1)

b <- 0.19

axis(1,c(0.5,max(xx)+0.5),line=0,labels=c("",""))
axis(2,c(0,0.1,0.2,0.3,0.4,0.5),line=0,labels=c(0,0.1,0.2,0.3,0.4,0.5))
#axis.break(2, b,style="slash")
#rect(0.5,b,max(xx)+0.5,b+0.005,border=FALSE,col='white')
text(xx[1,]+0.2,rep(-0.08,7), cex=1, offset=2,labels=ril5[,2],srt=55,xpd=TRUE)
legend('topright', bty='n', border='NA', lty=c(0,0),cex=1 ,c("Unique", "Sim_DHS_TSD", "Sim_DHS", "Sim_TSD", "Control"),fill=c("aquamarine3", "steelblue2" ,"sandybrown", "darkred", "dim gray"))



dev.off()

