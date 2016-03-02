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
python Excision_Number_In_RILs.py --input RIL230.sample.list

Sumamrize excision number in RILs
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#mPing   NumberOfEvents  NumberOfEvents_Footprint        NumberOfRILs    NumberOfRILs_Footprint  Footprint:RILs
#Chr10:11955070-11955072 1       0       1       0       Perfect:RIL185
#Chr10:13102744-13102746 0       0       0       0       0
#Chr10:16851282-16851284 3       2       4       2       16851260_D_20:RIL205;16851284_D_1-16851292_D_1-16851297_D_3-16851304_D_1:RIL237;Perfect:RIL92,RIL270
def read_excision_table(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'):
                unit = re.split(r'\t',line)
                if int(unit[3]) > 0 and not unit[5] == '0':
                    events = re.split(r';', unit[-1])
                    for e in events:
                        print e
                        footprint, ril = re.split(r':', e)
                        data[ril] += 1
    print '%s of RILs has excisions'  %(len(data.keys()))
    return data


#GN1
def read_rils(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                ril  = re.sub(r'GN', r'', unit[0])
                data[ril] =1
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

    rils = read_rils(args.input)
    excision_table = 'Excision_newpipe_version1.footprint.list.noPing.txt'
    excisions = read_excision_table(excision_table)
    ofile = open(re.sub(r'.txt', r'.rils.txt', excision_table), 'w')
    for ril in sorted(rils.keys(), key=int):
        ril = 'RIL%s' %(ril)
        print >> ofile, '%s\t%s' %(ril, excisions[ril])
    ofile.close()

if __name__ == '__main__':
    main()

