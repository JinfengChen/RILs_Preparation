plot_correlation_merged <- function(x, y, title, ...){
    plot(x, y, pch=18, col='orange', xlab= '', ylab= '', main= title, cex.main=1, xlim=c(0, 7), ylim=c(0,10))
    #regression line
    abline(lm(y~x), col="lightblue")
    reg <- lm(y~x)
    a <- round(reg$coefficients[[2]], 2)
    b <- round(reg$coefficients[[1]], 2)
    text(0, 8, pos=4, paste('Y = ', b, '+', paste(a, 'X', sep=''), sep=' '), cex=1)
    #corrlation
    cor <- cor.test(x, y)
    r2 <- round(cor$estimate[[1]], 2)
    p  <- cor$p.value
    if (p == 0){
        p = '2.2e-16'
    }else{
        p = signif(p, 2) 
    }
    text(0, 9, pos=4, paste('R2 =', r2, ', p-value <', p, sep=' '), cex=1)
    xpos <- 2.2
    ypos <- 5
    mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=0.8, col="black")
    mtext("copy number", side=1,font=1, at=xpos+2.1,line=3, cex=0.8, col="black")
    mtext("Excision number", side=2,font=1, at=ypos,line=3, cex=0.8, col="black")
}

plot_correlation_merged_parental_mping <- function(x, y, title, ...){
    plot(x, y, pch=18, col='orange', xlab= '', ylab= '', main= title, cex.main=1, xlim=c(120, 260), ylim=c(0, 10))
    #regression line
    abline(lm(y~x), col="lightblue")
    reg <- lm(y~x)
    a <- round(reg$coefficients[[2]], 2)
    b <- round(reg$coefficients[[1]], 2)
    text(120, 8, pos=4, paste('Y = ', b, '+', paste(a, 'X', sep=''), sep=' '), cex=1)
    #corrlation
    cor <- cor.test(x, y)
    r2 <- round(cor$estimate[[1]], 2)
    p  <- cor$p.value
    if (p == 0){
        p = '2.2e-16'
    }else{
        p = signif(p, 2) 
    }
    text(120, 9, pos=4, paste('R2 =', r2, ', p-value <', p, sep=' '), cex=1)
    xpos <- 180
    ypos <- 5
    mtext("Parental mPing copy number", side=1,font=1, at=xpos+10,line=3, cex=0.8, col="black")
    mtext("Excision number", side=2,font=1, at=ypos,line=3, cex=0.8, col="black")
    #mtext("number", side=2,font=1, at=ypos+183,line=3, cex=0.8, col="black") 
}


plot_correlation <- function(x, y, xlab, ylab, title, ...){
    plot(x, y, xlab= '', ylab= '', main= title)
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
    if (p == 0){
        p = '2.2e-16'
    }else{
        p = signif(p, 2) 
    }
    text(0, 130, pos=4, paste('R2 =', r2, ', p-value <', p, sep=' '), cex=1.4)
    xpos <- 2.6
    ypos <- 44
    mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1.4, col="black")
    mtext("copy number", side=1,font=1, at=xpos+1.6,line=3, cex=1.4, col="black")
    mtext("Unique homozygous", side=2,font=1, at=ypos,line=3, cex=1.4, col="black")
    mtext("mPing", side=2,font=3, at=ypos+42,line=3, cex=1.4, col="black")
    mtext("number", side=2,font=1, at=ypos+65,line=3, cex=1.4, col="black") 
}


pdf('ping_number_clean_excision_cor_plot.pdf')

x   <- read.table('RIL272.RIL_mPing_Ping_Excision.table.txt', header=TRUE)
x_n <- read.table('RIL272.RIL_mPing_Ping_Excision.narrow_range.table.txt', header=TRUE)

print("Parental mPing vs Excision, RIL272")
cor.test(x[,2], x[,11])
print("Ping vs Excision, RIL272")
cor.test(x[,9], x[,11])

print("Parental mPing vs Excision, RIL272_narrow_range")
cor.test(x_n[,2], x_n[,11])
print("Ping vs Excision, RIL272_narrow_range")
cor.test(x_n[,9], x_n[,11])

par(mfrow=c(2,2))
plot_correlation_merged_parental_mping(x[,2], x[,11], '272 RILs')
plot_correlation_merged(x[,9], x[,11], '272 RILs')
plot_correlation_merged_parental_mping(x_n[,2], x_n[,11], '123 narrow range mPing RILs')
plot_correlation_merged(x_n[,9], x_n[,11], '123 narrow range mPing RILs')

par(mfrow=c(1,1))
plot_correlation_merged_parental_mping(x[,2], x[,11], 'Parental mPing Copy Number', 'Excision Number', '272 RILs')
plot_correlation_merged(x[,9], x[,11], 'Ping Copy Number', 'Excision Number', '272 RILs')
plot_correlation_merged_parental_mping(x_n[,2], x_n[,11], 'Parental mPing Copy Number', 'Excision Number', '123 narrow range mPing RILs')
plot_correlation_merged(x_n[,9], x_n[,11], 'Ping Copy Number', 'Excision Number', '123 narrow range mPing RILs')

#parental mPing vs unique homozygous mPing given same number of ping
#default
#par(mfrow=c(3,2))
#for (ping in 2:7){
#     subx <- subset(x, Ping_Number==ping, select=c(Shared_HEG4, Unique_hom, Ping_Number))
#     title = paste(ping, 'Ping', sep=' ')
#     plot_correlation_merged_parental_mping(subx[,1], subx[,2], title)
#} 

#high quality
#par(mfrow=c(3,2))
#for (ping in 2:7){
#    subx <- subset(x_h, Ping_Number==ping, select=c(Shared_HEG4, Unique_hom, Ping_Number))
#    title = paste(ping, 'Ping', sep=' ')
#    plot_correlation_merged_parental_mping(subx[,1], subx[,2], title)
#}


dev.off()
