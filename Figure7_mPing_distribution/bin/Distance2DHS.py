#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
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



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def plot_R(prefix):
    cmd_R='''
error.bar <- function(x, y, upper, lower=upper, color,length=0.06,...){
     if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
     stop("vectors must be same length")
     arrows(x,y+upper, x, y-lower, col=color,angle=90, code=3, length=length, ...)
}

pdf("%s.pdf")

par(mar=c(6,4,4,2), cex=1.2)
dist <- read.table("%s.sum")
sim <- read.table("mping2DHS.simulation.sum")
plot(rev(dist[,4]), type='b', pch= 1,lwd = 2 , col="aquamarine3", xaxt='n', frame.plot = FALSE, ylim=c(0,0.3), ylab="Proportion", xlab="")
lines(rev(sim[,2]), type='b',pch= 20, cex=0.2,lwd = 2 , col="dim gray")
error.bar(1:length(sim[,2]), rev(sim[,2]), rev(sim[,3]), rev(sim[,3]), 'dim gray')
axis(1,seq(1:length(dist[,1])),line=0, labels=rep("",length(dist[,1])))
text(seq(1:length(dist[,1])),rep(-0.04,7), cex=1, offset=2,labels=rev(dist[,2]/1000),srt=0,xpd=TRUE)
mtext("Distance to DHS (kb)", side=1, cex=1.2, at=6, line=3)
legend('topright', bty='n', border='NA', lty= c(1,2), pch = c(1,20), cex=1 , lwd = 2 ,col=c("aquamarine3", "dim gray"), c("Unique", "Control"))
dev.off()
''' %(prefix, prefix)
    ofile = open ('%s.R' %(prefix), 'w')
    print >> ofile, cmd_R
    ofile.close()
    os.system('cat %s.R | R --slave' %(prefix))

def distance_sum_file(mping):
    data = defaultdict(lambda : int()) 
    inter = {
         -5:'1000',
         -4:'800',
         -3:'600',
         -2:'400',
         -1:'200',
         0:'0',
         1:'-200',
         2:'-400',
         3:'-600',
         4:'-800',
         5:'-1000',
    }
    tt = 0
    outer = 0
    for m in sorted(mping.keys()):
        l = mping[m]
        if 1 > 0:
            tt += 1
            if l > 1000 or l < -1000:
                continue
                outer += 1
            if l >= 800:
                data[-5] += 1
            elif l >= 600:
                data[-4] += 1
            elif l >= 400:
                data[-3] += 1
            elif l >= 200:
                data[-2] += 1
            elif l > 0:
                data[-1] += 1
            elif l == 0:
                data[0] += 1
            elif l >= -200:
                data[1] += 1
            elif l >= -400:
                data[2] += 1
            elif l >= -600:
                data[3] += 1
            elif l >= -800:
                data[4] += 1
            elif l >= -1000:
                data[5] += 1
    for r in sorted(data.keys(), key=int):
        print r, data[r]
        data[r] = float(data[r])/tt 
    return data

def distance_sum(mping, fname):
    data = defaultdict(lambda : int()) 
    inter = {
         -5:'1000',
         -4:'800',
         -3:'600',
         -2:'400',
         -1:'200',
         0:'0',
         1:'-200',
         2:'-400',
         3:'-600',
         4:'-800',
         5:'-1000',
    }
    tt = 0
    outer = 0
    for m in sorted(mping.keys()):
        l = mping[m]
        if 1 > 0:
            tt += 1
            if l > 1000 or l < -1000:
                continue
                outer += 1
            if l >= 800:
                data[-5] += 1
            elif l >= 600:
                data[-4] += 1
            elif l >= 400:
                data[-3] += 1
            elif l >= 200:
                data[-2] += 1
            elif l > 0:
                data[-1] += 1
            elif l == 0:
                data[0] += 1
            elif l >= -200:
                data[1] += 1
            elif l >= -400:
                data[2] += 1
            elif l >= -600:
                data[3] += 1
            elif l >= -800:
                data[4] += 1
            elif l >= -1000:
                data[5] += 1

    ofile = open(fname, 'w')
    for r in sorted(data.keys(), key=int):
        print >> ofile, '%s\t%s\t%s\t%s' %(r, inter[r], data[r], float(data[r])/int(tt))



def distance_table(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'): 
                unit = re.split(r'\t',line)
                mping= '%s_%s' %(unit[0], unit[3])
                distance = int(unit[18])
                if int(unit[12]) < int(unit[3]) and int(unit[13]) < int(unit[3]):
                    distance = 0 - distance
                data[mping] = distance
                #print infile, mping, distance
    return data

#-5      1000    73      0.0237862495927
#-4      800     72      0.0234604105572
def read_sum_file(sum_file):
    data = defaultdict(lambda : float())
    with open (sum_file, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t', line)
                data[int(unit[0])] = float(unit[3])
    return data


def simulation_sum_table(sum_files, fname):
    data = defaultdict(lambda : list())
    rank = [5,4,3,2,1,0,-1,-2,-3,-4,-5]
    total = 0
    for sum_file in sorted(sum_files):
        print sum_file
        sums = read_sum_file(sum_file)
        for r in rank:
            data[r].append(sums[r])
    ofile = open(fname, 'w')
    for r in sorted(data.keys(), key=int):
        print >> ofile, '%s\t%s\t%s' %(r, np.mean(data[r]), np.std(data[r])) 
    ofile.close()

#summary simulation for distance to DHS
def simulation_distance(simulation, dhs, sim_sum):
    gffs = glob.glob('%s/Simulate*.gff' %(simulation))
    sum_files = []
    for gff in sorted(gffs):
        bed = os.path.basename(re.sub(r'.gff', r'.closest.bed', gff))
        sums = os.path.basename(re.sub(r'.gff', r'.sum', gff))
        sum_files.append(sums)
        if not os.path.exists(bed):
            cmd = '/opt/bedtools/2.17.0-25-g7b42b3b/bin/bedtools closest -a %s -b %s -d > %s' %(gff, dhs, bed)
            os.system(cmd)
            mping2dhs = distance_table(bed)
            distance_sum(mping2dhs, sums)
    simulation_sum_table(sum_files, sim_sum)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    mping = '../input/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff'
    dhs   = '../input/GSM655033_Rice_Seedling_DHsites.MSU7.Corrected.gff'
    simulation = '../input/simulateV2_Random_TSD9mer_rilMat'
    cmd = '/opt/bedtools/2.17.0-25-g7b42b3b/bin/bedtools closest -a %s -b %s -d > mping2DHS.closest.bed' %(mping, dhs)
    if not os.path.exists('mping2DHS.closest.bed'):
        os.system(cmd)
    mping2dhs = distance_table('mping2DHS.closest.bed') 
    distance_sum(mping2dhs, 'mping2DHS.closest.sum')
    #simulation_distance(simulation, dhs, 'mping2DHS.simulation.sum')
 
    plot_R('mping2DHS.closest')
if __name__ == '__main__':
    main()

