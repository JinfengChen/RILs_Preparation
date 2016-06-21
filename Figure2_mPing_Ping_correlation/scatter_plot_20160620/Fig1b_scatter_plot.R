plot_correlation <- function(x, y, color='red', xpos=7, ypos=0.1, title="", ...){
  #regression line
  #abline(lm(y~x), col=color)
  #reg <- lm(y~x)
  #a <- round(reg$coefficients[[2]], 2)
  #b <- round(reg$coefficients[[1]], 2)
  #text(0, 120, pos=4, paste('Y = ', b, '+', paste(a, 'X', sep=''), sep=' '), cex=1)
  #corrlation
  cor <- cor.test(x, y)
  r2 <- round(cor$estimate[[1]], 2)
  p  <- cor$p.value
  if (p == 0){
    p = '2.2e-16'
  }else{
    p = signif(p, 2) 
  }
  text(xpos, ypos, pos=4, paste('R-squared =', r2, ', P =', p, sep=' '), cex=1.4, col=color) 
}

pdf("Fig1b_scatter_plot.pdf")
par(mar=c(5,5,4,2))
mping <- read.table("RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt", header=TRUE)
mping <- mping[ which(mping$Ping_Number<8),]
beeswarm(Unique_hom ~ Ping_Number, data = mping,
          pch = 18, col=c("cornflowerblue"),
          xlab = "", ylab = "",
          labels = c(0,1,2,3,4,5,6,7), axes=FALSE)
bxplot(Unique_hom ~ Ping_Number, probs=c(0.5), data = mping,col="darkgreen", add=TRUE)

plot_correlation(mping$Ping_Number, mping$Unique_hom, "cornflowerblue", 2, 120)


axis(2,seq(0, 140, by=20),line=0, labels=seq(0, 140, by=20), cex.axis=1.4)
axis(1,c(0.5, 8.5),line=0,labels=c("",""), cex.axis=1.4)
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