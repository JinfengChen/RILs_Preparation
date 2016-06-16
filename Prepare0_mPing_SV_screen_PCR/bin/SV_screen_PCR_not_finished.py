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
                if not data.has_key(unit[0]) and unit[1] != 'na':
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
                line = re.sub(r'RIL', r'', line)
                data.append(line) 
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    #../../Prepare0_mPing_SV_screen/bin/HEMP8_HEMP9.check.rils 
    rils_files = glob.glob('../../Prepare0_mPing_SV_screen/bin/*.check.rils')
    for f in rils_files:
        rils   = read_ril(f)
        result = './results/%s.PCR.results.txt' %(re.split(r'\.', os.path.split(f)[1])[0])
        if os.path.exists(result):
            header, rils_done = read_results(result)
            print os.path.split(f)[1]
            for ril in sorted(rils, key=int):
                if not rils_done.has_key('RIL%s' %(ril)):
                    print ril

if __name__ == '__main__':
    main()

