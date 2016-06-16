#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO
import scipy
from scipy import stats
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python NB_Ping_effect.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt --output Ping1_hom --ping A
    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


#Sample  Shared_HEG4     Shared_RILs     Shared  Unique  Unique_hom      Unique_het      Unique_som      Ping_Number     Ping_Code
#RIL1    159     50      209     28      19      9       0       5       DEFGH
#RIL2    113     21      134     12      10      2       0       2       CH
#RIL3    185     51      236     47      43      4       0       6       ABCEFG*
#RIL4    153     41      194     26      23      3       0       5       CDEFG
#RIL5    129     32      161     23      18      5       0       3       FGH
def NB_ping(infile, prefix, ping):
    data = defaultdict(lambda : defaultdict(lambda : list()))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'RIL') and not 'NA' in line: 
                unit = re.split(r'\t',line)
                #het = int(unit[6]) + int(unit[7])
                het = int(unit[5])
                #ping = 'A'
                if int(unit[8]) == 2 and ping in unit[9]:
                    data[2][1].append(het)
                    #print line
                    #print het
                elif int(unit[8]) == 2 and not ping in unit[9]:
                    data[2][0].append(het)
                    #print line
                    #print het
                elif int(unit[8]) == 3 and ping in unit[9]:
                    data[3][1].append(het) 
                elif int(unit[8]) == 3 and not ping in unit[9]:
                    data[3][0].append(het)
                elif int(unit[8]) == 4 and ping in unit[9]:
                    data[4][1].append(het)
                elif int(unit[8]) == 4 and not ping in unit[9]:
                    data[4][0].append(het)
                elif int(unit[8]) == 5 and ping in unit[9]:
                    data[5][1].append(het)
                elif int(unit[8]) == 5 and not ping in unit[9]:
                    data[5][0].append(het)
                elif int(unit[8]) == 6 and ping in unit[9]:
                    data[6][1].append(het)
                elif int(unit[8]) == 6 and not ping in unit[9]:
                    data[6][0].append(het)
                elif int(unit[8]) == 7 and ping in unit[9]:
                    data[7][1].append(het)
                elif int(unit[8]) == 7 and not ping in unit[9]:
                    data[7][0].append(het)
      
        ofile = open ('%s.table.txt' %(prefix), 'w')
        print >> ofile, 'Ping_Number\tCode\tHeterozygous_mPing\t'
        count = 0
        for ping_n in sorted(data.keys()):
            ping_w_nb_mean  = np.mean(data[ping_n][1])
            ping_w_nb_std   = np.std(data[ping_n][1])
            ping_wo_nb_mean = np.mean(data[ping_n][0])
            ping_wo_nb_std  = np.std(data[ping_n][0])
            for w in data[ping_n][1]:
                print >> ofile, '%s\t%s\t%s' %(ping_n, count, w)
            count += 1
            for wo in data[ping_n][0]:
                print >> ofile, '%s\t%s\t%s' %(ping_n, count, wo)
            count += 1
            w, pvalue = scipy.stats.mannwhitneyu(data[ping_n][1], data[ping_n][0])
            #w1, pvalue1 = scipy.stats.mannwhitneyu(data[ping_n][1], data[ping_n+1][1]) 
            #print '%s ping with NB_ping: n %s, mean %s, std %s' %(ping_n, len(data[ping_n][1]), ping_w_nb_mean, ping_w_nb_std)
            #print '%s ping without NB_ping: n %s, mean %s, std %s' %(ping_n, len(data[ping_n][0]), ping_wo_nb_mean, ping_wo_nb_std)
            #print 'Wilcoxon-Matt-Whitney test NB_ping vs Ping: %s, %s' %(w, pvalue)
            print '%s: %s/%s=%s' %(ping_n, ping_w_nb_mean, ping_wo_nb_mean, float(ping_w_nb_mean)/float(ping_wo_nb_mean))
            if ping_n < 7:
                w1, pvalue1 = scipy.stats.mannwhitneyu(data[ping_n][1], data[ping_n+1][1])
                print 'Wilcoxon-Matt-Whitney test low_ping vs high_ping: %s, %s' %(w1, pvalue1)
        ofile.close()

        Rcmd='''
pdf("%s.pdf")
par(mar=c(5,5,4,2))
library(beeswarm)
mping <- read.table("%s.table.txt", header=TRUE)
beeswarm(Heterozygous_mPing ~ Code, data = mping,
          pch = 16, col=c("cornflowerblue"),
          xlab = "", ylab = "",
          labels = c(2,3,4,5,6,7), axes=FALSE, ylim=c(0, 100), xlim=c(0.5, 12.5))
bxplot(Heterozygous_mPing ~ Code, data = mping, add=TRUE)
axis(2,seq(0, 100, by=20),line=0, labels=seq(0, 100, by=20), cex.axis=1.4)
axis(1,c(0.5, 12.5),line=0,labels=c("",""), cex.axis=1.4)
#x and y lab
xpos=5
ypos=30
#mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1.4, col="black")
#mtext("copy number", side=1,font=1, at=xpos+2.6,line=3, cex=1.4, col="black")
mtext("Unique homozygous", side=2,font=1, at=ypos,line=3, cex=1.4, col="black")
mtext("mPing", side=2,font=3, at=ypos+31,line=3, cex=1.4, col="black")
mtext("number", side=2,font=1, at=ypos+48,line=3, cex=1.4, col="black")
#x annotation
text(c(1,2,3,4,5,6,7,8,9,10,11,12), rep(-11, 12), offset=2,labels=c("+", "-", "+", "-","+", "-","+", "-","+", "-", "+", "-","+", "-","+", "-"),srt=0,xpd=TRUE, cex=1.4)
text(c(1.5,3.5,5.5,7.5,9.5,11.5), rep(-19, 12), offset=2,labels=c("2 Pings", "3 Pings", "4 Pings", "5 Pings", "6 Pings", "7 Pings"),srt=0,xpd=TRUE, cex=1.4)

dev.off()
''' %(prefix, prefix)
        ofile = open('%s.R' %(prefix), 'w')
        print >> ofile, Rcmd 
        ofile.close()
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-p', '--ping')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.ping:
        args.ping = 'H'

    NB_ping(args.input, args.output, args.ping)


if __name__ == '__main__':
    main()

