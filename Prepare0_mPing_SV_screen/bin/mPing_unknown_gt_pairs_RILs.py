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
python mPing_unknown_gt_pairs_RILs.py --input High_excision_mPing.pairs --matrix ../input/High_excision_csv_Ping

Find these RILs that have genotype not sure at pairs of high excision mPing. We only count these loci on HEG4 block.

Cases we treat as genotype known.
1. both mPing are present in the RILs (inversion between mPing will be not included in this case, as we can not see the difference between original and inversion using reads mapping).
2. both mPing are absence in the RILs (should be rare, but still exists. including these early events and co-excision)
3. one mPing excised and one mPing present (including these excision we valided).

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#../input/Bam.Core.blacklist
def read_blacklist(infile):
    data = defaultdict(lambda : str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]): 
                    data[unit[0]] = 1
    return data 

#Chr1.36267659 to Chr1_36267659_36267661
def id_reform(id1):
    unit = re.split(r'\.', id1)
    start = int(unit[1])
    end   = start + 2
    id2   = '%s_%s_%s' %(unit[0], start, end)
    return id2

#mPing1 mPing2  Distance        mPing1_primer   mPing2_primer   Note
#Chr1.36267659   Chr1.36270511   2852    HEMP1   HEMP2
#Chr1.6806761    Chr1.6816415    9654    HEMP3   HEMP4
def read_mping_pairs(infile):
    data = defaultdict(lambda : list())
    pair = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                unit = re.split(r'\t',line)
                pairs= '%s_%s' %(unit[3], unit[4])
                pair.append(pairs)
                id1 = id_reform(unit[0])
                id2 = id_reform(unit[1])                
                data[pairs] = [id1, id2, unit[2]]
    return data, pair


#Chr1:29494572-29494575  133,134,158,185,271,87
def read_mping_table(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = 1
    return data

#Chr5:25474861-25474863,Genotype_Bin,Pseudo_mPing_status_up,Pseudo_mPing_status_down,Ref_mPing_status,Excision_Code,Ping_Number,Ping_Code
#RIL1,HEG4,covered,covered,unknown,Insertion,5,DEFGH
#RIL2,HEG4,unknown,covered,unknown,Insertion,2,CH
def read_mping_gt(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'RIL'): 
                unit = re.split(r',',line)
                ril = re.sub(r'RIL', r'', unit[0])
                data[ril] = unit
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-m', '--matrix')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
   
    if not args.output:
        args.output = 'High_excision_mPing'

    ofile_sum = open('%s.summary' %(args.output), 'w')
    print >> ofile_sum, 'Pairs\tTotal\tNB\tHEG4\tHEG4_Not_Done'
    blacklist = read_blacklist('../input/Bam.Core.blacklist.230')
    pairs, pair_list = read_mping_pairs(args.input)
    for pair in pair_list:
        mping1 = pairs[pair][0]
        mping2 = pairs[pair][1]
        gt1     = '%s/%s.matrix.csv' %(args.matrix, mping1)
        gt2     = '%s/%s.matrix.csv' %(args.matrix, mping2)
        mping_gt1 = read_mping_gt(gt1)
        mping_gt2 = read_mping_gt(gt2)
        stat = [0,0,0,0,0] # total, nb, heg4, not done
        ofile = open('%s.check.list' %(pair), 'w')
        print >> ofile, '%s,Genotype_Bin,Pseudo_mPing_status_up,Pseudo_mPing_status_down,Ref_mPing_status,Excision_Code,Ping_Number,Ping_Code,%s,Genotype_Bin,Pseudo_mPing_status_up,Pseudo_mPing_status_down,Ref_mPing_status,Excision_Code,Ping_Number,Ping_Code,mPing1_status,mPing2_status,Validation,Blastlist' %(mping1, mping2)
        for ril in sorted(mping_gt1.keys(), key=int):
            black = 'Blackout' if blacklist.has_key('RIL%s' %(ril)) else 'Okay'
            stat[0] += 1
            if mping_gt1[ril][1] == mping_gt2[ril][1] and mping_gt1[ril][1] == 'HEG4':
                stat[2] += 1
                mping1_ins, mping2_ins = [0, 0]
                mping1_exc, mping2_exc = [0, 0]
                if mping_gt1[ril][2] == mping_gt1[ril][3] and mping_gt1[ril][2] == 'covered' and mping_gt1[ril][4] == 'clipped':
                    mping1_ins = 1
                if mping_gt2[ril][2] == mping_gt2[ril][3] and mping_gt2[ril][2] == 'covered' and mping_gt2[ril][4] == 'clipped':
                    mping2_ins = 1
                if mping_gt1[ril][2] == mping_gt1[ril][3] and mping_gt1[ril][2] == 'clipped' and mping_gt1[ril][4] == 'covered':
                    mping1_exc = 1
                if mping_gt2[ril][2] == mping_gt2[ril][3] and mping_gt2[ril][2] == 'clipped' and mping_gt2[ril][4] == 'covered':
                    mping2_exc = 1
                if mping1_ins == 1 and mping2_ins == 1:
                    print >> ofile, '%s,%s,Insertion,Insertion,finished' %(','.join(mping_gt1[ril]), ','.join(mping_gt2[ril]))
                elif mping1_exc == 1 and mping2_exc == 1:
                    print >> ofile, '%s,%s,Excision,Excision,finished' %(','.join(mping_gt1[ril]), ','.join(mping_gt2[ril]))
                elif mping1_ins == 1 and mping2_exc == 1:
                    print >> ofile, '%s,%s,Insertion,Excision,finished' %(','.join(mping_gt1[ril]), ','.join(mping_gt2[ril]))
                elif mping1_exc == 1 and mping2_ins == 1:
                    print >> ofile, '%s,%s,Excision,Insertion,finished' %(','.join(mping_gt1[ril]), ','.join(mping_gt2[ril]))
                else:
                    stat[3] += 1
                    status1 = 'NA'
                    if mping1_ins == 1:
                        status1 = 'Insertion'
                    elif mping1_exc == 1:
                        status1 = 'Excision'
                    status2 = 'NA'
                    if mping2_ins == 1:
                        status2 = 'Insertion'
                    elif mping2_exc == 1:
                        status2 = 'Excision'
                    print >> ofile, '%s,%s,%s,%s,not_done,%s' %(','.join(mping_gt1[ril]), ','.join(mping_gt2[ril]), status1, status2, black)
            elif mping_gt1[ril][1] == mping_gt2[ril][1] and mping_gt1[ril][1] == 'NB':
                print >> ofile, '%s,%s,NA,NA,NA,%s' %(','.join(mping_gt1[ril]), ','.join(mping_gt2[ril]), black)
                stat[1] += 1
            else:
                print >> ofile, '%s,%s,NA,NA,NA,%s' %(','.join(mping_gt1[ril]), ','.join(mping_gt2[ril]), black)
                stat[4] += 1
        ofile.close()
        print >> ofile_sum, '%s\t%s' %(pair, '\t'.join(map(str, stat)))
    ofile_sum.close()
    
if __name__ == '__main__':
    main()

