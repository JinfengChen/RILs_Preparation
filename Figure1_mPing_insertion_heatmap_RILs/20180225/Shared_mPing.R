plot_heatmap <- function(ind){
    header <- colnames(ind)
    ind_m <- as.matrix(ind)
    ind_m_melt <- melt(ind_m)
    base_size = 10
    p <- ggplot(data=ind_m_melt, aes(x=X2, y=X1)) + geom_tile(aes(fill=value), colour="snow2") + scale_fill_gradient(low='lightblue', high="steelblue") + theme_grey(base_size = base_size) + labs(x = "",y = "") + theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank()) + xlim(header)  + ylim(row.names(ind))
    return(p) 
}

pdf("Shared_mPing.pdf", width=7, height=5)
par(mar=c(4,5,10,4))

library(reshape)
library(ggplot2)
library(Rmisc)
ind_p <- read.table("Shared_Parental_0.3.matrix", header=TRUE)
ind_s <- read.table("Shared_RIL_0.1.matrix", header=TRUE)
ind_p_p <- plot_heatmap(ind_p)
ind_s_p <- plot_heatmap(ind_s)

library("ggpubr")
ggarrange(ind_p_p, ind_s_p, ncol=2, widths = c(2, 9))

dev.off()
