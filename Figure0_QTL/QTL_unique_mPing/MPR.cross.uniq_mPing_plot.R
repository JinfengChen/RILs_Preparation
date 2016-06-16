library(qtl)
pdf ("MPR.cross.uniq_mPing_plot.pdf", width=10, height=4)
par(mar=c(5,5,4,2))
## step0. read and write the data of cross
if (1){
   read.cross("qtlcart",dir="./",file="MPR.cross.uniq.cro",mapfile="MPR.cross.uniq.map") -> cross
}else{
   read.cross("csv",dir="./",file="MPR.cross.uniq") -> cross
}

## step1. do single QTL analysis
## 1.1 marker regression
QTL.mr <- scanone(cross,pheno.col=7,method="mr")

## 1.2.1 estimate the threshold of LOD using permutation test
#operm <- scanone(cross,pheno.col=7,n.perm=1000, verbose=FALSE)
#LOD <- summary(operm,alpha=0.05) ## get LOD threshold
#QTL.mr.test <- summary(QTL.mr,perms=operm,alpha=0.05,format="allpeak",pvalues=TRUE)
chr <- 1:12
plot(QTL.mr, chr=c(chr), lodcolumn=1, main='', ylab="LOD score", col="black", lwd=1)
abline(h=3.30, lty=2, col="orange")
#iplotMScanone(QTL.mr)
#ping_pos <- c("17.924119", "23.439904", "136.703516", "125.669806", "33.193178", "48.885324", "73.890293", "87.340104")
ping_pos <- c(17.924119, 23.439904, 366.76+25*2+136.703516, 1027.97+25*6+125.669806, 1300.73+25*8+33.193178, 1300.73+25*8+48.885324, 1300.73+25*8+73.890293, 892.32+25*5+87.340104)
ping_lod <- c(6.81, 4.32, 0.17, 0.65, 3.51, 4.74, 1.04, 3.63e-6)
ping_name <- c("Ping_A", "Ping_B", "Ping_C", "Ping_D", "Ping_E", "Ping_F", "Ping_G", "Ping_H")
for (i in 1:8){
    segments(ping_pos[i], 7.2, ping_pos[i], ping_lod[i], col='gray', lty=2, xpd=TRUE)
    text(ping_pos[i], 7.4, labels=ping_name[i], adj=0, srt=45, xpd=TRUE)
}
dev.off()



