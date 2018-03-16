pdf("ping_pong_high_strains.heatmap.pdf")
library(reshape)
library(ggplot2)
library(Rmisc)
header <- c("NB", "HEG4", "A119", "A123", "IRIS_313.15904", "B160", "B235", "OtherStrains")
ind <- read.table("ping_pong_high_strains.matrix", row.names=1, header=TRUE)
ind_m <- as.matrix(ind)
ind_m_melt <- melt(ind_m)
base_size = 10
ggplot(data=ind_m_melt, aes(x=X2, y=X1)) + geom_tile(aes(fill=value), colour="snow2") + scale_fill_gradient(low='lightblue', high="steelblue") + theme_grey(base_size = base_size) + labs(x = "",y = "") + scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "none",axis.ticks = element_blank(), axis.text.x = element_text(size = base_size *0.8, angle = 40, hjust = 1, vjust=1, colour = "grey50")) + xlim(header)  + ylim(row.names(ind))

dev.off()
