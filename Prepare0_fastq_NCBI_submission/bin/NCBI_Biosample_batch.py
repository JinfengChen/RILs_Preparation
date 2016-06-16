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

#RIL1    1       /rhome/cjinfeng/Rice/RIL/Illumina_correct/RIL1_0/RIL1_0_CGTACG_FC153L5_p1.fq.gz
def read_list(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                ril  = re.sub(r'RIL', r'', unit[0])
                data[ril] = 1
    return data

def NCBI_biosample(rils, tsv):
    ofile = open(tsv, 'a')
    for ril in sorted(rils.keys(), key=int):
        #print ril
        print >> ofile, 'RIL%s\tnot collected\tPRJNA316308\tOryza Sativa\tnot collected\tNipponbare x HEG4 RIL%s\tTemperate Japonica\tnot collected\tSeedling\tnot collected\tLeaf\tSusan Wessler\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected' %(ril, ril) 
    ofile.close()

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

    #cp Plant.1.0.tsv RILs_272.NCBI_biosample.tsv
    os.system('cp ../input/Plant.1.0.tsv RILs_272.NCBI_biosample.tsv')
    rils = read_list(args.input) 
    NCBI_biosample(rils, 'RILs_272.NCBI_biosample.tsv')

if __name__ == '__main__':
    main()

