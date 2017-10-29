#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob

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

def runjob(script, lines):
    #cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 30 --lines %s --interval 120 --task 12 --mem 100G --time 20:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-r', '--ref')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
   
    if not args.ref:
        #args.ref = '/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/input/bwa_0.7/MSU_r7.Pseudo_mPing_RILs.fa'
        #args.ref = '/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Pseudogenome/MSU_r7.Pseudo_mPing.fa'
        #args.ref = '/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Pseudogenome/MSU_r7.Pseudo_mPing_RILs.fa'
        args.ref = '/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Pseudogenome/MSU_r7.Pseudo_mPing_415.fa'
    fqs = glob.glob('%s/RIL*/*_1.fq' %(args.input))
    ofile = open('RIL_bwa.sh', 'w')
    for fq1 in sorted(fqs):
        fq1 = os.path.abspath(fq1)
        fq2 = re.sub(r'_1', r'_2', fq1)
        fq1_dirs = re.split(r'/', fq1)
        ril = fq1_dirs[-2]
        cmd = 'perl /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Identification/bin/step1_Mapping.pl -ref %s -1 %s -2 %s -min 0 --max 500 -cpu 12 --tool bwa --project %s --verbose' %(args.ref, fq1, fq2, os.path.abspath(ril))
        if not os.path.exists('{}.bam'.format(ril)):
            print >> ofile, cmd
    ofile.close()
    
    runjob('RIL_bwa.sh', 1)
if __name__ == '__main__':
    main()

