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
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0 and len(args.frequency)
    except:
        usage()
        sys.exit(2)

    sum_dict = defaultdict(lambda : defaultdict(lambda : int()))
    mpings = gff_parser(args.gff)
    frq    = readtable(args.frequency)
    for mping in sorted(mpings.keys()):
        mping_idx = '%s:%s-%s' %(mpings[mping][0][0], mpings[mping][0][3], mpings[mping][0][4])
        mping_tp  = mpings[mping][1]['type']
        tp_idx    = mping_tp[:3]
        if frq.has_key(mping_idx):
            #print '%s\t%s\t%s' %(mping_idx, mping_tp, frq[mping_idx])
            if int(frq[mping_idx]) < 10:
                sum_dict['2.Less than 10 RILs'][tp_idx] += 1 
            else:
                sum_dict['1.More than 10 RILs'][tp_idx] += 1
        else:
            #print '%s\t%s\t%s' %(mping_idx, mping_tp, 'NA')
            sum_dict['3.Absence in RILs'][tp_idx] += 1
    for t in sorted(sum_dict.keys()):
        temp = [t]
        for h in ['hom', 'het', 'som']:
            temp.append(sum_dict[t][h])
        print '\t'.join([str(i) for i in temp])

if __name__ == '__main__':
    main()

