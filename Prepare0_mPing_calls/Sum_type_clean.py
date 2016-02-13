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
python Sum_type_clean.py --prefix RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean

Summary types of mPing call: unique, parental or shared in RILs
    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

def write_file(fn, lines):
    ofile = open(fn, 'w')
    print >> ofile, lines
    ofile.close()


def parse_gff(infile):
    num = 0
    if not os.path.isfile(infile) or not os.path.getsize(infile):
        return 0 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                num += 1 
    return num


def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def parse_table(infile):
    parental = 0
    ril      = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'): 
                unit = re.split(r'\t',line)
                #print unit[-1]
                if unit[-1] == 'Parental':
                    parental += 1
                elif unit[-1] == 'RIL':
                    ril      += 1
    return [parental, ril]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--prefix')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.prefix) > 0
    except:
        usage()
        sys.exit(2)

    prefix = args.prefix
    shared_mping_file = '%s.shared_mping.ril.frequency' %(prefix)
    unique_mping_gff  = '%s.unique_mPing.gff' %(prefix)
    type_file = '%s.type.summary' %(prefix) 
    unique_mping     = parse_gff(unique_mping_gff)
    parental, shared = parse_table(shared_mping_file)  
    #print unique_mping, parental, shared
    write_file(type_file, 'Unique\t%s\nShared\t%s\nParental\t%s' %(unique_mping, shared, parental)) 
if __name__ == '__main__':
    main()

