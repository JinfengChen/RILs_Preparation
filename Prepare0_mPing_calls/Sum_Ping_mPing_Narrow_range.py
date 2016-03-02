#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python Unique_mPing.py --input RIL275_RelocaTEi.CombinedGFF.characterized.gff

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Pings	Ping_Code	RIL
#0	NA	RIL129
#1	C	RIL230
def summary_unique(ping_code, unique, ping_number_sum, ping_single_sum):
    data = defaultdict(lambda : list())
    #ril_mping =  read_mping_table(mping_table)
    ril_unique=  read_unique_table(unique)
    with open (ping_code, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Pings'):
                unit   = re.split(r'\t',line)
                ril_id = re.sub(r'RIL', 'RIL', unit[2])
                data[ril_id] = [unit[0], unit[1]]
    sum_ping_number_hom = defaultdict(lambda : list())
    sum_ping_number_som = defaultdict(lambda : list())
    sum_ping_number_m = defaultdict(lambda : list())
    sum_ping_single_hom = defaultdict(lambda : list())
    sum_ping_single_som = defaultdict(lambda : list())
    sum_ping_single_m = defaultdict(lambda : list())
    for ril in ril_unique.keys():
        if not data[ril][0] == 'NA':
            sum_ping_number_hom[data[ril][0]].append(ril_unique[ril][0])
            sum_ping_number_som[data[ril][0]].append(ril_unique[ril][1])
            sum_ping_number_m[data[ril][0]].append('%s:%s:%s' %(ril, ril_unique[ril][0], ril_unique[ril][1]))
        if len(data[ril][1]) == 1:
            sum_ping_single_hom[data[ril][1]].append(ril_unique[ril][0])
            sum_ping_single_som[data[ril][1]].append(ril_unique[ril][1])
            sum_ping_single_m[data[ril][1]].append('%s:%s:%s' %(ril, ril_unique[ril][0], ril_unique[ril][1]))
    ofile0 = open(ping_number_sum, 'w')
    print >> ofile0, '#Ping\tUnique_mPing_hom_mean\tUnique_mPing_hom_std\tUnique_mPing_het_mean\tUnique_mPing_het_std\tUnique_mPing_hom_list\tUnique_mPing_het_list\tUnique_mPing_RILs'
    ofile1 = open(ping_single_sum, 'w')
    print >> ofile1, 'Ping\tUnique_mPing_hom_mean\tUnique_mPing_hom_std\tUnique_mPing_het_mean\tUnique_mPing_het_std\tUnique_mPing_hom_list\tUnique_mPing_het_list\tUnique_mPing_RILs'
    for n_ping in sorted(sum_ping_number_hom.keys(), key = int):
        values_hom = map(int, sum_ping_number_hom[n_ping])
        values_som = map(int, sum_ping_number_som[n_ping])
        #print 'n_ping', values_hom, values_som
        print >> ofile0, n_ping, np.mean(values_hom), np.std(values_hom), np.mean(values_som), np.std(values_som), len(values_hom), ','.join(map(str, values_hom)), ','.join(map(str, values_som)), ','.join(sum_ping_number_m[n_ping]) 
    for s_ping in sorted(sum_ping_single_hom.keys()):
        values_hom = map(int, sum_ping_single_hom[s_ping])
        values_som = map(int, sum_ping_single_som[s_ping])
        #print 's_ping', values_hom, values_som
        print >> ofile1, s_ping, np.mean(values_hom), np.std(values_hom), np.mean(values_som), np.std(values_som), len(values_hom), ','.join(map(str, values_hom)), ','.join(map(str, values_som)), ','.join(sum_ping_single_m[n_ping])
    ofile0.close()
    ofile1.close()
    return data

#Sample  InsertSize      Size_STD        Depth   Mapped_Depth    Mapped_Rate     Non_Ref_mPing   Confident       Evidence_From_One_End   Characterized   Homozygous      Heterozygous    Somatic
#RIL1    141.39  25.74   6.35    6.19    0.98    336     228     108     226     209     17      0
def read_mping_table(infile):
    data  = defaultdict(lambda : list())
    r     = re.compile(r'RIL(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Sample'):
                unit  = re.split(r'\t', line)
                #print unit[10], unit[11], unit[12]
                data[unit[0]] = [unit[10], int(unit[11]) + int(unit[12])]
                #print data[unit[0]]
    return data

#Sample  Shared_HEG4     Shared_RILs     Shared  Unique  Unique_hom      Unique_het      Unique_som
#RIL1    151     45      196     30      20      10      0
def read_unique_table(infile):
    data  = defaultdict(lambda : list())
    r     = re.compile(r'RIL(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Sample'):
                unit  = re.split(r'\t', line)
                #print unit[10], unit[11], unit[12]
                if int(unit[1]) >= 173 and int(unit[1]) <= 212:
                    data[unit[0]] = [unit[5], str(int(unit[6]) + int(unit[7]))]
                #print data[unit[0]]
    return data


def main():
    parser = argparse.ArgumentParser()
    #parser.add_argument('-t', '--table')
    parser.add_argument('-c', '--code')
    parser.add_argument('-o', '--output')

    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    #if not args.table:
    #    args.table     = '../RIL275_RelocaTEi.summary_clean.table'
    if not args.code:
        args.code      = 'RIL275_RelocaTE.sofia.ping_code.table'

    #unique = 'RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.txt'
    if not args.output:
        args.output = 'RIL275_RelocaTEi.CombinedGFF.characterized.clean'
    prefix = 'RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean'
    unique = '%s.mping.shared_unique_table.txt' %(prefix)
    ping_number_sum = '%s.ping_number.summary' %(args.output)
    ping_single_sum = '%s.ping_single.summary' %(args.output)
    summary_unique(args.code, unique, ping_number_sum, ping_single_sum)

if __name__ == '__main__':
    main()

