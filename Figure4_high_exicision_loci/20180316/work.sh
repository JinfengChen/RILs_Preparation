echo "identification of high exicison mPing using mPing excision frequency"
cat mping.excision_events.distr.R | R --slave
cat mping.excision_events.binomial_test.R | R --slave

