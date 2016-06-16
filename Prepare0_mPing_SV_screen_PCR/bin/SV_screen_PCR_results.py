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
python SV_screen_PCR_results.py


    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def read_results(infile):
    header = ''
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'RIL'): 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] =[unit[1], unit[2]]
            elif line.startswith(r'HEMP'):
                header = line
    return header, data


def read_ril(infile):
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                data.append(line) 
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
 
    rils    = read_ril('RILs.list')
    results = glob.glob('./*.PCR.results.txt')
    for f in results:
        header, data = read_results(f)
        new_file = re.sub(r'.txt', '.rils.txt', f)
        ofile = open(new_file, 'w')
        print >> ofile, header
        for ril in rils:
            if data.has_key(ril):
                print >> ofile, '%s,%s,%s' %(ril, data[ril][0], data[ril][1])
            else:
                print >> ofile, '%s,--,--' %(ril)
        ofile.close()

if __name__ == '__main__':
    main()

