#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import HTSeq
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)


def plot_TSD_profile_R(profile, prefix):
    cmd_R='''
pdf("%s.pdf")
par(mar=c(6,4,4,2), cex=1.2)
dist <- read.table('%s')
dist <- subset(dist, V1<=2500 & V1>=500)
dist_sim <- read.table('%s.simulation.profile.sum')
dist_sim <- subset(dist_sim, V1<=2500 & V1>=500)
plot(rev(dist[,3]/1000000), type='l', pch= 1,lwd = 2 , col="aquamarine3", xaxt='n', frame.plot = FALSE, ylim=c(10000/1000000,25000/1000000), ylab="Normalized DNase reads", xlab="")
lines(rev((dist_sim[,2]*10/1000000)-0.005), type='b',pch= 20, cex=0.2,lwd = 2 , col="dim gray")
#error.bar(1:length(dist[,2]), rev(dist[,2]), rev(dist[,3]), rev(dist[,3]), 'dim gray')
axis(1,seq(0, 2000, by=200),line=0, labels=rep("",11))
text(seq(0, 2000, by=200), rep(0.008, 11), cex=1, offset=2,labels=seq(-1, 1, by=0.2), srt=0,xpd=TRUE)
mtext("Distance to TSD (kb)", side=1, cex=1.2, at=1000, line=3)
legend('topright', bty='n', border='NA', lty= c(1,2), pch = c(1,20), cex=1 , lwd = 2 ,col=c("aquamarine3", "dim gray"), c("Unique", "Control"))
dev.off()
''' %(prefix, profile, prefix)
    ofile = open('%s.R' %(prefix), 'w')
    print >> ofile, cmd_R
    ofile.close()
    os.system('cat %s | R --slave' %('%s.R' %(prefix)))

def profile_TSD(gff, bam, halfwinwidth, fragmentsize, total, tsd_profile_sum):
    ##gff
    #Chr1    RelocaTE        mPing   10046640        10046642        .       .       .       Strains=HEG4_2;GT=heterozygous
    #Chr1    RelocaTE        mPing   10508960        10508962        .       .       .       Strains=RIL;GT=heterozygous
    tsdpos = []
    with open (gff, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith('Chr'): 
                unit = re.split(r'\t',line)
                #print unit[0], unit[3], unit[4]
                tsdpos.append([unit[0], unit[3], unit[4]])

    sortedbamfile = HTSeq.BAM_Reader(bam)
    ofile = open(tsd_profile_sum, 'w')
    profile_ALL = np.zeros( 2*halfwinwidth, dtype='i' )
    for p in tsdpos:
        #define empty vector of profile
        #profile = numpy.zeros( 2*halfwinwidth, dtype='i' )
        #define 1000 bp window around TSD
        window = HTSeq.GenomicInterval( p[0], int(p[1]) - halfwinwidth, int(p[1]) + halfwinwidth, "." )
        #for all the aligned reads in window
        for almnt in sortedbamfile[ window ]:
            #set reads to fragmentsize, which might be 200 or 300
            almnt.iv.length = fragmentsize
            #convert reads coordinate on chromosome to coordinate in window
            start_in_window = almnt.iv.start - int(p[1]) + halfwinwidth
            end_in_window   = almnt.iv.end   - int(p[1]) + halfwinwidth
            #make sure coordinate within the window range becaue we extended reads to 200 bp
            start_in_window = max( start_in_window, 0 )
            end_in_window = min( end_in_window, 2*halfwinwidth )
            if start_in_window >= 2*halfwinwidth or end_in_window < 0:
                continue
            #add depth to profile
            #profile[start_in_window : end_in_window] += 1
            profile_ALL[start_in_window : end_in_window] += 1
       
        #mping = '%s_%s' %(p[0], p[1])
        #depth = '\t'.join(map(str, profile))
        #print >> ofile, '%s\t%s' %(mping, depth)

    for i in range(len(profile_ALL)):
        print >> ofile, '%s\t%s\t%s' %(i, profile_ALL[i], float(profile_ALL[i])/total)
    ofile.close()


def read_profile_file(infile):
    data = defaultdict(lambda : float())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = float(unit[2])
    return data

def simulation_profile_table(profile_files, fname):
    data = defaultdict(lambda : list())
    for profile_file in sorted(profile_files):
        if not os.path.exists(profile_file):
            continue
        profiles = read_profile_file(profile_file)
        for r in sorted(profiles.keys()):
            data[r].append(profiles[r])
    ofile = open(fname, 'w')
    for r in sorted(data.keys(), key=int):
        print >> ofile, '%s\t%s\t%s' %(r, np.mean(data[r]), np.std(data[r])) 
    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam')
    parser.add_argument('-g', '--gff')
    parser.add_argument('-s', '--simulation')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0
    except:
        usage()
        sys.exit(2)

    if not args.bam:
        args.bam = '../input/DHS.unique.bam'
    if not args.output:
        args.output = 'mping_TSD'
    if not args.simulation:
        args.simulation = '../input/simulateV2_Random_TSD9mer_rilMat'

    halfwinwidth = 1500
    fragmentsize = 200
    #total = 60745783.00/1000000 ## nucleosome
    total = 23299296/1000000 #DHS unique
    gff = args.gff
    bam = args.bam
    tsd_profile = '%s.profile' %(args.output)
    tsd_profile_sum = '%s.profile.sum' %(args.output)    

    #unique mPing
    print 'Calculating unique mPing DHS profile...'
    if not os.path.exists(tsd_profile_sum):
        print 'profiling %s' %(gff)
        profile_TSD(gff, bam, halfwinwidth, fragmentsize, total, tsd_profile_sum)
        #plot_TSD_profile_R(tsd_profile_sum, args.output)
    else:
        print '%s exists' %(tsd_profile_sum)
    print 'unique mPing DHS profile done'    

    
    #simulation mPing
    print 'Calculating simulation mPing DHS profile...'
    sim_gffs = glob.glob('%s/Sim*.gff' %(args.simulation))
    sim_profile_sums = []
    for sim_gff in sorted(sim_gffs):
        print sim_gff
        sim_profile_sum = os.path.basename(re.sub(r'.gff', r'.profile.sum', sim_gff))
        if not os.path.exists(sim_profile_sum):
            print 'profiling %s' %(gff)
            try:
                profile_TSD(sim_gff, bam, halfwinwidth, fragmentsize, total, sim_profile_sum)
            except:
                continue
        else:
            print '%s exists' %(sim_profile_sum)
        sim_profile_sums.append(sim_profile_sum)
    print 'Simulation mPing DHS profile done' 

    
    #summary simulation data
    print 'Summary and plot profile'
    sim_profile_table = '%s.simulation.profile.sum' %(args.output)
    simulation_profile_table(sim_profile_sums, sim_profile_table)
    plot_TSD_profile_R(tsd_profile_sum, args.output)
    print 'ALL Done'

if __name__ == '__main__':
    main()

