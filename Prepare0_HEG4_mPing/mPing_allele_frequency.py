#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser

def usage():
    test="name"
    message='''
python mPing_allele_frequency.py --frequency ../input/mping.ril.frquency --gff ../input/HEG4.ALL.mping.non-ref.gff
python mPing_allele_frequency.py --frequency ../input/mping.ril.frquency --gff ../input/HEG4.ALL.mping.non-ref.gff --subgff HEG4.ALL.mping.non-ref.AF0.1.gff > HEG4.ALL.mping.non-ref.allele.frq

Generate allele frequency table of mPing from gff

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr9    7356518 7356520 Chr9:7356518_7356520    +       1       0.00363636363636
def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                index= '%s:%s-%s' %(unit[0], unit[1], unit[2])
                data[index] = unit[5]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gff')
    parser.add_argument('-f', '--frequency')
    parser.add_argument('--subgff')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0 and len(args.frequency)
    except:
        usage()
        sys.exit(2)

    mpings = gff_parser(args.gff)
    frq    = readtable(args.frequency)
    ofile  = ''
    if args.subgff:
        ofile = open(args.subgff, 'w')
    for mping in sorted(mpings.keys()):
        mping_idx = '%s:%s-%s' %(mpings[mping][0][0], mpings[mping][0][3], mpings[mping][0][4])
        if frq.has_key(mping_idx):
            print '%s\t%s' %(mping_idx, frq[mping_idx])
            if args.subgff:
                if float(frq[mping_idx]) >= 0:
                    print >> ofile, '\t'.join(mpings[mping][0])
        else:
            print '%s\t%s' %(mping_idx, 'NA')
    if args.subgff:
        ofile.close()

if __name__ == '__main__':
    main()

