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


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#mPing   NumberOfEvents  NumberOfEvents_Footprint        NumberOfRILs    NumberOfRILs_Footprint  Footprint:RILs
#Chr1:29494572-29494574  6       5       11      7       29494564_D_17:RIL87;29494572_D_
def read_high_excision(infile):
    data = defaultdict(lambda : str())
    r = re.compile(r'(\w+\:\d+)\-\d+')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if r.search(unit[0]):
                    mping = r.search(unit[0]).groups(0)[0]
                    data[mping] = unit[1]
    return data

#LOC_Os01g01010.1        TBC domain containing protein,expressed
def read_msu_anno(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = unit[1]
    return data

#Chr1_1715119    Chr1    1715117 1715119 Five_primer_UTR LOC_Os01g03980.1        -2170
def read_nearby_gene(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = '%s:%s' %(unit[1], unit[2])
                data[mping] = [unit[5], unit[6]]
    return data

#mPing   Excision        Distance
#Chr1.10903901   0       322234
def read_distance(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith(r'Chr'): 
                unit = re.split(r'\t',line)
                mping = re.sub(r'\.', r':', unit[0])
                #print mping, unit[2]
                data[mping] = unit[2]
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

    high_excisions   = read_high_excision('Excision_newpipe_version1.footprint.high.txt')
    HEG4_distance    = read_distance('Excision_distance_HEG4.matrix_events.1.txt')
    RIL_distance     = read_distance('Excision_distance_RIL.matrix_events.1.txt')

    print 'High excision mPing\tDistance to closest mPing loci in HEG4\tDistance to closest mPing loci in RILs'
    for mping in sorted(high_excisions.keys()):
        print '%s\t%d kb\t%d kb' %(mping, float(HEG4_distance[mping])/1000, float(RIL_distance[mping])/1000)

if __name__ == '__main__':
    main()

