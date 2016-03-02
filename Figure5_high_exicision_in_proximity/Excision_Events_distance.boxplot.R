pvalue <- function(x1, y1, x2, y2, top, p){
  star <- '*'
  if (p > 0.05) {star <- 'n.s.'}
  if (p < 0.001){ star <- '**'}
  if (p < 0.0001){ star <- '***'}
  segments(x1,y1,x1,top)
  segments(x1,top,x2,top)
  segments(x2,top,x2,y2)
  #segments(x1-0.2,y1,x1+0.2,y1)
  #segments(x2-0.2,y2,x2+0.2,y2)
  xt <- min(x1,x2)+abs(x2-x1)/2
  yt <- top*1.04
  p <- paste("P-value = ", p, sep="")
  text(xt,yt,p, cex=1, xpd=TRUE)
} 

pdf("Excision_Events_distance.boxplot.pdf", height=7, width=4)
par(mar=c(7,4,4,2))
high <- read.table("Excision_Events_distance.boxplot.high.txt")
all <- read.table("Excision_Events_distance.boxplot.all.txt")
#all <- t(all)
#high <- t(high)
total <- length(all[,1]) + length(high[,1])
data <- c(all[,1], high[,1])
mat <- matrix(c(data, 0), ncol=2)
df <- data.frame(values = mat[1:total], vars = rep(c("all","high"), times = c(length(all[,1]),length(high[,1]))))
values = mat[1:total]
vars = rep(c("all","high"), times = c(length(all[,1]),length(high[,1])))
boxplot(values ~ vars, data = df, ylim=c(0, 1200000), axes=FALSE, ylab="Distance (kb)")
axis(1,c(0.6, 2.4),line=0,labels=c("",""))
#mtext("All Parental mPing")
text(1, -240000, "Parental mPing", xpd=TRUE, srt=35)
text(2, -240000, "High frequency excision mPing", xpd=TRUE, srt=35)
axis(2,seq(0, 1200000, by=200000), line=0, labels=seq(0, 1200, by=200))
wilcox.test(all[,1], high[,1])
pvalue(1, 1230000, 2, 400000, 1250000, 1.40e-05)
dev.off()
