plot_correlation <- function(x, y, color='red', xpos=7, ypos=0.1, title="", ...){
  #regression line
  abline(lm(y~x), col=color)
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
  text(xpos, ypos, pos=4, paste(title, ",", 'P =', p, sep=' '), cex=1.4, col=color) 
}

pdf("Fig1d_scatter_plot_ping_rt.pdf")
par(mar=c(5,5,4,2))

library(beeswarm)
ping_rt <- read.table("Ping_transcription.txt", header=TRUE, sep="\t")

beeswarm(ORF1 ~ Ping_Number, data = ping_rt,
pch = 17, col=c("blue"),
xlab = "", ylab = "",
labels = c(1,2,3,4,5,6,7), ylim=c(0, 0.35), axes=FALSE)

plot_correlation(ping_rt$Ping_Number, ping_rt$ORF1, "blue", 3.5, 0.015, "ORF1/Actin")

par(new=TRUE)

beeswarm(Tpase ~ Ping_Number, data = ping_rt,
pch = 15, col=c("red"),
xlab = "", ylab = "",
labels = c(1,2,3,4,5,6,7), ylim=c(0, 0.35), axes=FALSE)

plot_correlation(ping_rt$Ping_Number, ping_rt$Tpase, "red", 2, 0.27, "Tpase/Actin")

#x and y lab
axis(2,seq(0, 0.35, by=0.05),line=0, labels=seq(0, 0.35, by=0.05), cex.axis=1.4)
axis(1,c(0.5, 7.5),line=0,labels=c("",""), cex.axis=1.4)
xpos=3
ypos=0.2
mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1.4, col="black")
mtext("copy number", side=1,font=1, at=xpos+1.65,line=3, cex=1.4, col="black")
mtext("Relative expression level", side=2,font=1, at=ypos,line=3, cex=1.4, col="black")
#x annotation
text(c(1,2,3,4,5,6,7), rep(-0.03, 7), offset=2,labels=c(1,2,3,4,5,6,7),srt=0,xpd=TRUE, cex=1.4)

dev.off()