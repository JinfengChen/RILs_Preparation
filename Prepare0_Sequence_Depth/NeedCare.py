#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import os.path
import commands
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python inf_scripts/NeedCare.py > NeedCare.bam.inf

Generate summary for sequence coverage and percentage of NA for SNPs. 
RIL	Qual	Coverage	PercentOfNA
1	solexa/illumina	6.59613	0.347223633525
2	solexa/illumina	5.53489	0.411498988481
3	solexa/illumina	6.87857	0.327584454934
4	solexa/illumina	6.05613	0.320720903809

PercentOfNA:
Not Done is not analyzed using QTL pipeline
    '''
    print message
#
#GN252_? 38282338        101     3866516138      10.3939
#GN34_?  36416222        101     3678038422      9.8872

#Sample  #Read   Depth   Mapped_Depth    Mapped_rate     Dupli_rate      Insert_median   Map_quality     GC_percent      Coverage_mapped Coverage(1-5X)  BamFile
#1       23383042        6.24616643385   6.09326612105   0.975520935215  0.0652200000325 132     49.93    39.99  91.98   91.98%;83.92%;72.33%;60.48%;48.21%      RIL1_0_CGTACG_FC153L5.recal.bam
def fq_cvg(infile):
    data = defaultdict(float)
    s = re.compile(r'GN(\d+)\_\?')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if not line.startswith('Samp;e') and len(line) > 2: 
                unit = re.split(r'\t',line)
                #m = s.search(unit[0])
                #ril = m.groups(0)[0] if m else '000'
                ril = unit[0]
                data[ril] = unit[2]
    return data            
    
'''
guess fastq format
'''
def fq_qual(rils):
    data = defaultdict(str)
    for ril in rils.keys():
        fq = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_fastq_correct/GN' + str(ril) + '_1.fq'
        if os.path.isfile(fq):
            cmd = 'perl /rhome/cjinfeng/BigData/software/bin/fastqformatdetect.pl ' + fq
            status, fqformat = commands.getstatusoutput(cmd)
            unit = re.split(r' ', fqformat) 
            data[ril]=unit[3]
        else:
            data[ril]='sanger' 
    return data
 
'''
GN1     GN10    GN100   GN101   GN102   GN103   GN104   GN105   GN106   GN107   GN108   GN109   GN11    GN110
0100021547A     NA      NA      NA      A       A       A       NA      NA      NA      NA      A       A
0100031071A     NA      G       A       A       A       A       G       G       G       G       NA      A
'''
def snpmatrix(infile):
    data = defaultdict(lambda : int())
    data1 = defaultdict(lambda: float())
    rils = []
    total = 0
    s = re.compile(r'GN(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            m = s.search(line)
            if line.startswith('GN') and m: 
                unit = re.split(r'\t',line)
                rils = [a.replace('GN','') for a in unit]
                print 'Number of RILs:', len(rils)
            elif len(line) > 2:
                total = total + 1
                unit = re.split(r'\t',line)
                unit = unit[1:]
                #print 'Number of genotype:', len(unit)
                count = 0
                for snp in unit: 
                    count = count + 1
                    #if int(rils[count-1]) == 40:
                    #    print 'RIL40 SNP:', snp
                    if snp == 'NA':
                        #if int(rils[count-1]) == 40:
                        #    print 'Called NA:', snp
                        #    print data[rils[count-1]]
                        data[rils[count-1]] = int(data[rils[count-1]]) + 1
    
    for ril in data.keys():
        narate = float(data[ril])/float(total)
        data1[ril] = narate
        #print ril, data[ril], total, narate
    return data1



'''
Sample  Heading Days    Plant Height (cm) in Field      Biomass Number of Tillers       Single Plant (Grain Yield) (g)
GN-1    103     92      120.2   22      42.0
'''
def trait(infile):
    data = defaultdict(str)
    s = re.compile(r'GN-(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            m = s.search(line)
            if line.startswith('GN') and m: 
                unit = re.split(r'\t',line)
                ril = m.groups(0)[0]
                if not data.has_key(ril):
                    data[ril] = 1
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    if args.input is None:
        args.input = '1'

    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
   
    rils = trait('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/trait/May28_2013.RIL.trait.table.QTL.trait.txt')
    snps = snpmatrix('./NB.RILs.dbSNP.SNPs.RILs')
    qual = fq_qual(rils) 
    #cvg  = fq_cvg('./RIL.fastq.fastq.stat')   
    #cvg  = fq_cvg('./inf_list/RIL.bam.unique.stat')
    cvg = fq_cvg('RILs_ALL_bam_correct_merged.sorted.summary')
    print "RIL\tQual\tCoverage\tPercentOfNA" 
    for ril in sorted(rils.keys(), key=int):
        coverage = cvg[ril] if cvg[ril] else 'NA'
        quality  = qual[ril] if qual[ril] else 'NA'
        if snps.has_key(ril):
            if snps[ril] > 0:
                print "%s\t%s\t%s\t%s" %(ril, quality, coverage, snps[ril])
        else:
            if qual.has_key(ril):
                print "%s\t%s\t%s\t%s" %(ril, quality, coverage, 'Not Done')
            else:
                print "%s\t%s\t%s\t%s" %(ril, quality, coverage, 'No sequence')

if __name__ == '__main__':
    main()

