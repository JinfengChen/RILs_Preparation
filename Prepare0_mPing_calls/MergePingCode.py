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
python MergeRILinfo.py --input info.list
Merge information from several list together using RIL ID.

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


#Sample  Shared_HEG4     Shared_RILs     Shared  Unique  Unique_hom      Unique_het      Unique_som
#RIL1    151     45      196     30      20      10      0
def merge_file(infile, ping_code, prefix):
    data  = defaultdict(lambda : list())
    r     = re.compile(r'RIL(\d+)')
    ofile = open(prefix, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Sample'):
                unit  = re.split(r'\t', line)
                ril   = re.sub(r'RIL', r'RIL', unit[0])
                if ping_code.has_key(ril):
                    unit.extend(ping_code[ril])
                else:
                    unit.extend(['NA','NA'])
                print >> ofile, '\t'.join(unit)
            else:
                print >> ofile, '%s\tPing_Number\tPing_Code' %(line)
    ofile.close()
    return data




def read_ping_code(infile):
    data  = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Pings'):
                unit   = re.split(r'\t',line)
                ril_id = re.sub(r'RIL', 'RIL', unit[2])
                data[ril_id] = [unit[0], unit[1]]
    return data

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

    prefix = '%s.ping_code.txt' %(os.path.splitext(args.input)[0])
    print prefix
    ping_code = read_ping_code('RIL275_RelocaTE.sofia.ping_code.table')
    merge_file(args.input, ping_code, prefix)

if __name__ == '__main__':
    main()

