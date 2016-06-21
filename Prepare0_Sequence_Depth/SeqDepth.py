#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python SeqDepth.py --input ../input/RIL.bam.stat

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#1	sanger	6.34862487634	0.28627553457
#2	sanger	5.33498941935	0.352280202204
def readinf(infile):
    data = defaultdict(lambda: float())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = unit[3]
    return data

#Sample  #Read   Average Total   Depth   Mapped_Depth    Mapped_rate     #Library        FileName
#GN1_?   23383054        101     2361688454      6.34862487634   6.19321574194   0.975520819479  1       ../../input/fastq/Bam/RIL1_0_CGTACG_FC153L5.recal.bam

#Sample  #Read   Depth   Mapped_Depth    Mapped_rate     Dupli_rate      Insert_median   Map_quality     GC_percent      Coverage_mapped Coverage(1-5X)  BamFile
#1       23383042        6.24616643385   6.09326612105   0.975520935215  0.0652200000325 132     49.93    39.99  91.98   91.98%;83.92%;72.33%;60.48%;48.21%      RIL1_0_CGTACG_FC153L5.recal.bam
def readtable(infile, inf):
    depth = defaultdict(lambda: float())
    mapped = defaultdict(lambda: float())
    s = re.compile(r'GN(\d+)\_')
    ofile = open('Depth.table', 'w')
    print >> ofile, 'RILs\tDepth\tMapped_Depth\tNA%'
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Sample'): 
                unit = re.split(r'\t',line)
                #m = s.search(unit[0])
                #ril = m.groups(0)[0] if m else 0
                ril = unit[0]
                depth[ril] += float(unit[2])
                mapped[ril] += float(unit[3])
    for r in sorted(depth.keys(), key=int):
        print >> ofile, '%s\t%s\t%s\t%s' %(r, depth[r], mapped[r], inf[r])
    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    inf = readinf('RILs_ALL_272.needcare.txt')
    readtable(args.input, inf)

if __name__ == '__main__':
    main()

