pdf("ping_effect_hom_ratio.pdf", height=5, width=6)
library(reshape)
library(ggplot2)
library(Rmisc)
header <- c("2 Pings","3 Pings","4 Pings","5 Pings","6 Pings","7 Pings")
ind <- read.table("Ping_effect_hom_ratio.txt", row.names=1, header=TRUE)
ind_m <- as.matrix(ind)
ind_m_melt <- melt(ind_m)
ind_m_melt[,3] <- log(ind_m_melt[,3])
base_size = 10
#ggplot(data=ind_m_melt, aes(x=X2, y=X1)) + geom_tile(aes(fill=value), colour="snow2") + scale_fill_gradient(low='lightblue', high="steelblue") + theme_grey(base_size = base_size) + labs(x = "",y = "") + scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "none",axis.ticks = element_blank(), axis.text.x = element_text(size = base_size *0.8, angle = 40, hjust = 1, vjust=1, colour = "grey50")) + xlim(header)  + ylim(row.names(ind))
ggplot(data=ind_m_melt, aes(x=X2, y=X1)) + geom_tile(aes(fill=value), colour="snow2") + scale_fill_gradient(low="steelblue", high="red") + theme_grey(base_size = base_size) + labs(x = "",y = "") + scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0))

dev.off()
