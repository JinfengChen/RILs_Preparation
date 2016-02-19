#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
Chr1    RIL10_0 RelocaTE        1132975 1132977 .       .       .       ID=mPing_1;Strain=RIL10_0;TSD=TAA;
'''
def readtable(infile):
    data = defaultdict(lambda: defaultdict(lambda: int))
    rils = defaultdict(int)
    inf  = defaultdict()
    p = re.compile(r'Strain=(\w+);')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = '%s:%s' %(unit[0], unit[3])
                m = p.search(unit[8])
                strain = m.groups(0)[0] if m else 'NA'
                data[mping][strain] = 1
                rils[strain] =1
                inf[mping] = [unit[0], unit[3], unit[4]]
    total = len(rils.keys())
    #print total
    for m in data.keys():
        count = len(data[m].keys())
        frq = float(count)/total
        print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(inf[m][0], inf[m][1], inf[m][2], m, '+', str(count), str(frq))
    


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

