#P=909/555/138=0.012 we expected only 50% of the RILs (138=275/2) have these 555 mPing from HEG4.
#P=341/419/115=0.007 we exptect only 50% of the RILS (115=130/2) have these 419 mPing from HEG4
#binom.test(5, 136, p=0.012)
#pvalue=0.036 when we observed 5 or more independent events for one locus

#update 20150628
#p=909/(555*(275/2)*2)=0.0059. 275/2 means only 50% of RILs have these 555 mPing from HEG4. *2 means diploid, which have two mping per loci.
#pvalue=0.001 when we observed 5 or more independent events for one locus
#binom.test(5, 137, p=0.0059)
#pvalue=0.009 when we observed 5 or more independent events for one locus
#binom.test(4, 137, p=0.0059)


#update 20150830
#P=341/(419*(230/2)*2)=0.0035 we exptect only 50% of the RILS (115=230/2) have these 419 mPing from HEG4
#p-value=0.0003039
#p-value=5.855e-05, 20151028
#binom.test(5, 115, p=0.0035)
#p-value=0.002781
#p-value=0.0007,    20151028
#binom.test(4, 115, p=0.0035)

#update 20160229
#P=342/(413*(272/2)*2)=0.0030
#p-value=6.311e-05
#binom.test(5, 136, p=0.0030)
#p-value=0.0008
#binom.test(4, 136, p=0.0030)

#update 20180316
#https://en.wikipedia.org/wiki/Binomial_test#Large_samples
#P=344/(466*(272/2)*2)=0.0027 #diploid
#P=344/(466*(272/2))=0.0054   #haploid
#p-value=0.0009201 or 9.2e-04
binom.test(5, 136, p=0.0054)
#p-value=0.006593  or 6.6e-03
binom.test(4, 136, p=0.0054)
binom.test(3, 136, p=0.0054)
binom.test(2, 136, p=0.0054)
binom.test(1, 136, p=0.0054)

#p-value=3.85e-05
binom.test(5, 136, p=0.0027)
#p-value=0.0005456 or 5.4e-04
binom.test(4, 136, p=0.0027)
binom.test(3, 136, p=0.0027)
binom.test(2, 136, p=0.0027)
binom.test(1, 136, p=0.0027)

#the rate is excision per locus per plant.
#the test is for the given rate, what the the chance of one mPing can have more than 5 excision events, which also mean in 5 rils.
