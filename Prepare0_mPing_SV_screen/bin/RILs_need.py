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

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#RIL1,NB,unknown,clipped,unknown,NB_GT,5,DEFGH,RIL1,NB,clipped,unknown,unknown,NB_GT,5,DEFGH,NA,NA,NA,Okay
def readtable(infile):
    ofile = open(re.sub(r'\.list', '.rils', infile), 'w')
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'RIL'): 
                unit = re.split(r',',line)
                if not unit[18] == 'finished' and not unit[1] == 'NB' and not unit[9] == 'NB' and not unit[-1] == 'Blackout':
                    #print line
                    ril = re.sub(r'RIL', r'', unit[0])
                    data[ril] = 1
                    print >> ofile, 'RIL%s' %(ril)
    ofile.close()
    return data

#mping.excision.RILs.list.uniq
def read_list(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                ril = re.sub(r'RIL', r'', unit[0])
                data[ril] = 1
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    #try:
    #    len(args.input) > 0
    #except:
    #    usage()
    #    sys.exit(2)
 
    #read list of RIL already have DNA for high excision validation
    dna_0 = read_list('mping.excision.RILs.list.uniq')

    #read matrix and extract RIL that need DNA do SV screen
    all_pairs = glob.glob('./*.check.list')
    dna_1 = defaultdict(lambda : str())
    for f in all_pairs:
        dna_1.update(readtable(f))

    #output RIL we need DNA
    for ril in sorted(dna_1.keys(), key=int):
        if not dna_0.has_key(ril):
            print 'RIL%s' %(ril)
    

if __name__ == '__main__':
    main()

