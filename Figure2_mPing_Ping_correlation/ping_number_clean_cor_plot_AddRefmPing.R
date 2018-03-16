plot_correlation_merged <- function(x, y, xlab, ylab, title, ...){
    plot(x, y, xlab= '', ylab= '', main= title, cex.main=1, xlim=c(0, 7), ylim=c(0, 140), pch=18)
    #regression line
    abline(lm(y~x), col="red")
    reg <- lm(y~x)
    a <- round(reg$coefficients[[2]], 2)
    b <- round(reg$coefficients[[1]], 2)
    text(0, 120, pos=4, paste('Y = ', b, '+', paste(a, 'X', sep=''), sep=' '), cex=1)
    #corrlation
    cor <- cor.test(x, y)
    print(cor)
    r2 <- round(cor$estimate[[1]], 2)
    p  <- cor$p.value
    print(p)
    #if (p == 0){
    #    p = '2.2e-16'
    #}else{
    p = signif(p, 2) 
    print(p)
    #}
    text(0, 130, pos=4, paste('r =', r2, ', P =', p, sep=' '), cex=1)
    xpos <- 2.2
    ypos <- 44
    mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=0.8, col="black")
    mtext("copy number", side=1,font=1, at=xpos+2.1,line=3, cex=0.8, col="black")
    mtext("Unique homozygous", side=2,font=1, at=ypos,line=3, cex=0.8, col="black")
    mtext("mPing", side=2,font=3, at=ypos+66,line=3, cex=0.8, col="black")
    mtext("number", side=2,font=1, at=ypos+103,line=3, cex=0.8, col="black") 
}

plot_correlation_merged_parental_mping <- function(x, y, title, ...){
    plot(x, y, pch=18, col='orange', xlab= '', ylab= '', main= title, cex.main=1, xlim=c(0, 300), ylim=c(0, 140))
    #regression line
    abline(lm(y~x), col="lightblue")
    reg <- lm(y~x)
    a <- round(reg$coefficients[[2]], 2)
    b <- round(reg$coefficients[[1]], 2)
    text(10, 100, pos=4, paste('Y = ', b, '+', paste(a, 'X', sep=''), sep=' '), cex=1)
    #corrlation
    cor <- cor.test(x, y)
    r2 <- round(cor$estimate[[1]], 2)
    p  <- cor$p.value
    #if (p == 0){
    #    p = '2.2e-16'
    #}else{
        p = signif(p, 2) 
    #}
    text(10, 120, pos=4, paste('r =', r2, ', P =', p, sep=' '), cex=1)
    xpos <- 60
    ypos <- 44
    mtext("Parental", side=1,font=1, at=xpos+10,line=3, cex=0.8, col="black")
    mtext("mPing", side=1,font=3, at=xpos+70,line=3, cex=0.8, col="black")
    mtext("copy number", side=1,font=1, at=xpos+148,line=3, cex=0.8, col="black")
    mtext("Unique homozygous", side=2,font=1, at=ypos,line=3, cex=0.8, col="black")
    mtext("mPing", side=2,font=3, at=ypos+118,line=3, cex=0.8, col="black")
    #mtext("number", side=2,font=1, at=ypos+183,line=3, cex=0.8, col="black") 
}


plot_correlation <- function(x, y, xlab, ylab, title, ...){
    plot(x, y, xlab= '', ylab= '', main= title, xlim=c(0, 7), ylim=c(0, 140))
    #regression line
    abline(lm(y~x), col="red")
    reg <- lm(y~x)
    a <- round(reg$coefficients[[2]], 2)
    b <- round(reg$coefficients[[1]], 2)
    text(0, 120, pos=4, paste('Y = ', b, '+', paste(a, 'X', sep=''), sep=' '), cex=1.4)
    #corrlation
    cor <- cor.test(x, y)
    r2 <- round(cor$estimate[[1]], 2)
    p  <- cor$p.value
    #if (p == 0){
    #    p = '2.2e-16'
    #}else{
        p = signif(p, 2) 
    #}
    text(0, 130, pos=4, paste('r =', r2, ', P =', p, sep=' '), cex=1.4)
    xpos <- 2.6
    ypos <- 44
    mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1.4, col="black")
    mtext("copy number", side=1,font=1, at=xpos+1.6,line=3, cex=1.4, col="black")
    mtext("Unique homozygous", side=2,font=1, at=ypos,line=3, cex=1.4, col="black")
    mtext("mPing", side=2,font=3, at=ypos+42,line=3, cex=1.4, col="black")
    mtext("number", side=2,font=1, at=ypos+65,line=3, cex=1.4, col="black") 
}


pdf('ping_number_clean_cor_plot_AddRefmPing.pdf')
x <- read.table("RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.Ref_mPing.txt", header=T)
x <- x[ which(x$Ping_Number<8),]
x_n <- read.table("RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.Ref_mPing.narrow_range.txt", header=T)
x_h <- read.table("RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.Ref_mPing.high_depth.txt", header=T)
x_hn <- read.table("RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.Ref_mPing.high_narrow.txt", header=T)
hist(x[,2], xlab = 'HEG4 mPing copy number', main='')
plot(density(x[,2]), xlim=c(0, 300), xpd=TRUE, xlab = 'Parental mPing copy number', main='')
abline(v=c(199, 243), col=c('blue', 'blue'))
hom <- x[,6]
hom_n <- x_n[,6] 
hom_h <- x_h[,6]
hom_hn <- x_hn[,6]
#som <- x[,7]+x[,8]
ping <- x[,9]
ping_n <- x_n[,9]
ping_h <- x_h[,9]
ping_hn <- x_hn[,9]
#cor.test(ping, hom)
#cor.test(ping_n, hom_n)
#cor.test(ping_h, hom_h)
#cor.test(ping_hn, hom_hn)

par(mfrow=c(2,2))
plot_correlation_merged(ping, hom, 'Ping Copy Number', 'Unique homozygous mPing number', 'Default')
plot_correlation_merged(ping_n, hom_n, 'Ping Copy Number', 'Unique homozygous mPing number', 'Parental mPing range (199-243)')
plot_correlation_merged(ping_h, hom_h, 'Ping Copy Number', 'Unique homozygous mPing number', 'High quality (depth>=6X and mapped>=90%)')
plot_correlation_merged(ping_hn, hom_hn, 'Ping Copy Number', 'Unique homozygous mPing number', 'High quality and Parental mPing range')

par(mfrow=c(1,1))
plot_correlation(ping, hom, 'Ping Copy Number', 'Unique homozygous mPing number', 'Default')
plot_correlation(ping_n, hom_n, 'Ping Copy Number', 'Unique homozygous mPing number', 'Parental mPing range (199-243)')
plot_correlation(ping_h, hom_h, 'Ping Copy Number', 'Unique homozygous mPing number', 'High quality (depth>=6X and mapped>=90%)')
plot_correlation(ping_hn, hom_hn, 'Ping Copy Number', 'Unique homozygous mPing number', 'High quality and Parental mPing range')

#parental mPing vs unique homozygous mPing given same number of ping
#default
par(mfrow=c(3,2))
for (ping in 2:7){
     subx <- subset(x, Ping_Number==ping, select=c(Shared_HEG4, Unique_hom, Ping_Number))
     title = paste(ping, 'Pings', sep=' ')
     plot_correlation_merged_parental_mping(subx[,1], subx[,2], title)
} 

#high quality
par(mfrow=c(3,2))
for (ping in 2:7){
    subx <- subset(x_h, Ping_Number==ping, select=c(Shared_HEG4, Unique_hom, Ping_Number))
    title = paste(ping, 'Pings', sep=' ')
    plot_correlation_merged_parental_mping(subx[,1], subx[,2], title)
}


dev.off()
