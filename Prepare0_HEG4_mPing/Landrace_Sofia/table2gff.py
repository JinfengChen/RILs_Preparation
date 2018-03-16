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
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines 2 --interval 120 --task 1 --mem 15G --time 100:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#EG4_2	mping	TTA	Chr1.10903903	47	0	homozygous
#Chr1	EG4_2	transposable_element_attribute	2927	2927	.	.	.	ID=Chr1.2927.spanners;avg_flankers=43.5;spanners=0;type=homozygous;TE=mping;TSD=TTA
def readtable(infile):
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                chrs, start = re.split(r'\.', unit[3])
                temp = []
                temp.append(chrs)
                temp.append(unit[0])
                temp.append('transposable_element_attribute')
                temp.append(start)
                temp.append(start) 
                temp.append('.')
                temp.append('.')
                temp.append('.')
                temp.append('ID={}.spanners;avg_flankers={};spanners={};type={};TE={};TSD={}'.format(unit[3], unit[4], unit[5], unit[6], unit[1], unit[2]))
                print '\t'.join(temp)

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
    
    readtable(args.input)

if __name__ == '__main__':
    main()

