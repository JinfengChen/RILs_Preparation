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

Generate number of unique mPing list for each RILs
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
def summary_unique(unique_mping_d, ping_code, ping_number_sum, ping_single_sum):
    data = defaultdict(lambda : list())
    with open (ping_code, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Pings'):
                unit   = re.split(r'\t',line)
                ril_id = re.sub(r'RIL', '', unit[2])
                data[ril_id] = [unit[0], unit[1]]
    sum_ping_number = defaultdict(lambda : list())
    sum_ping_number_m = defaultdict(lambda : list())
    sum_ping_single = defaultdict(lambda : list())
    sum_ping_single_m = defaultdict(lambda : list())
    for ril in unique_mping_d.keys():
        if not data[ril][0] == 'NA':
            sum_ping_number[data[ril][0]].append(unique_mping_d[ril])
            sum_ping_number_m[data[ril][0]].append('RIL%s:%s' %(ril, unique_mping_d[ril]))
        if len(data[ril][1]) == 1:
            sum_ping_single[data[ril][1]].append(unique_mping_d[ril])
            sum_ping_single_m[data[ril][1]].append('RIL%s:%s' %(ril, unique_mping_d[ril]))
    ofile0 = open(ping_number_sum, 'w')
    print >> ofile0, 'Class\t#Ping\t#Unique_mPing_mean\t#Unique_mPing_std\tSample_Size\tValues'
    ofile1 = open(ping_single_sum, 'w')
    print >> ofile1, 'Class\t#Ping\t#Unique_mPing_mean\t#Unique_mPing_std\tSample_Size\tValues'
    for n_ping in sorted(sum_ping_number.keys(), key = int):
        values = map(int, sum_ping_number[n_ping])
        #ril_mping = ','.join(sum_ping_number_m[n_ping])
        print >> ofile0, n_ping, np.mean(values), np.std(values), len(values), ','.join(map(str, sum_ping_number[n_ping])), ','.join(sum_ping_number_m[n_ping])

    for s_ping in sorted(sum_ping_single.keys()):
        values = map(int, sum_ping_single[s_ping])
        #ril_mping = ','.join(sum_ping_single_m[s_ping])
        print >> ofile1, s_ping, np.mean(values), np.std(values), len(values), ','.join(map(str, sum_ping_single[s_ping])), ','.join(sum_ping_single_m[s_ping])
    ofile0.close()
    ofile1.close()
    return data

#Chr1    RIL231_0        transposable_element_attribute  4228091 4228092 +       .       .       ID=Chr1.4228092.spanners;Strain=RIL231_0;
def unique_mping(infile, overlap_ref_d, overlap_ril_d, output):
    data  = defaultdict(lambda : int())
    r     = re.compile(r'RIL(\d+)_\d+')
    ofile = open(output, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                unit  = re.split(r'\t', line)
                index = '%s:%s_%s_%s' %(unit[1], unit[0], unit[3], unit[4])
                if not overlap_ref_d.has_key(index) and not overlap_ril_d.has_key(index):
                    print >> ofile, line
                    ril_id = r.search(unit[1]).groups(0)[0] if r.search(unit[1]) else 'NA'
                    data[ril_id] += 1
    ofile.close()
    return data

#Chr1    RIL231_0        transposable_element_attribute  4228091 4228092 +       .       .       ID=Chr1.4228092.spanners;Strain=RIL231_0;
def parse_overlap(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                unit = re.split(r'\t', line)
                if not unit[1] == unit[10] and unit[3] == unit[12] and unit[4] == unit[13]:
                    index       = '%s:%s_%s_%s' %(unit[1], unit[0], unit[3], unit[4])
                    data[index] = 1
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-r', '--reference')
    parser.add_argument('-w', '--window')
    parser.add_argument('-c', '--code')
    parser.add_argument('-o', '--output')

    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.reference:
        args.reference = 'HEG4.ALL.mping.non-ref.gff'
    if not args.window:
        args.window    = 0
    if not args.output:
        args.output    = '%s.unique_mPing.gff' %(os.path.splitext(args.input)[0])
    if not args.code:
        args.code      = 'RIL275_RelocaTE.sofia.ping_code.table'


    prefix = os.path.splitext(os.path.splitext(args.input)[0])[0]   
    

    overlap_ref     = '%s.overlap_ref_NB' %(os.path.splitext(args.input)[0])
    bed_overlap_ref = 'bedtools intersect -wb -f 1 -a %s -b %s > %s' %(args.input, args.reference, overlap_ref)
    if not os.path.exists(overlap_ref):
        os.system(bed_overlap_ref)
    overlap_ril     = '%s.overlap_ref_ril' %(os.path.splitext(args.input)[0])
    bed_overlap_ril = 'bedtools intersect -wb -f 1 -a %s -b %s > %s' %(args.input, args.input, overlap_ril)
    if not os.path.exists(overlap_ril):
        os.system(bed_overlap_ril)
    #overlap_ref_d   = parse_overlap(overlap_ref)
    #overlap_ril_d   = parse_overlap(overlap_ril)
    #unique_mping_d  = unique_mping(args.input, overlap_ref_d, overlap_ril_d, args.output)
    #ping_number_sum = '%s.ping_number.unique.summary' %(os.path.splitext(args.input)[0])
    #ping_single_sum = '%s.ping_single.unique.summary' %(os.path.splitext(args.input)[0])
    #summary_unique(unique_mping_d, args.code, ping_number_sum, ping_single_sum)

if __name__ == '__main__':
    main()

