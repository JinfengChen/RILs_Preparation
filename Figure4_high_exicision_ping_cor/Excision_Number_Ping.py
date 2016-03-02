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
python Excision_Number_Ping.py --excision Excision_newpipe_version1.footprint.list.noPing.rils.txt --ping RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.narrow_range.txt

Read list of mPing/Ping number and excision numbers. Combine them together in output file.

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


#Sample  Shared_HEG4     Shared_RILs     Shared  Unique  Unique_hom      Unique_het      Unique_som      Ping_Number     Ping_Code
#RIL1    150     45      195     30      20      10      0       5       DEFGH
#or
#RIL1    3
def read_table(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'RIL'):
                unit = re.split(r'\t',line)
                ril  = re.sub(r'RIL', r'', unit[0])
                data[ril] = unit
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--excision')
    parser.add_argument('-p', '--ping')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.excision) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = 'RIL230.RIL_mPing_Ping_Excision.narrow_range.table.txt'

    mping_table    = read_table(args.ping)
    excision_table = read_table(args.excision)
    ofile = open(args.output, 'w')
    print >> ofile, 'Sample\tShared_HEG4\tShared_RILs\tShared\tUnique\tUnique_hom\tUnique_het\tUnique_som\tPing_Number\tPing_Code\tExcision'
    for ril in sorted(mping_table.keys()):
        print >> ofile, '%s\t%s' %('\t'.join(mping_table[ril]), excision_table[ril][1])
    ofile.close()

if __name__ == '__main__':
    main()

