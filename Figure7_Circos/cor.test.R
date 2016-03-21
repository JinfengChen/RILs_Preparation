ril <- read.table('GFF.RIL.histogram.txt')
het <- read.table('GFF.Somatic.histogram.txt')
cor.test(ril[,3], het[,4])

