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
python Sum_Shared_mPing_Ref.py --gff RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.Shared.gff

Summary Reference_only and shared mPing for each RIL
 
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr1    not.give        transposable_element_attribute  4220010 4220012
def read_gff(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'): 
                unit = re.split(r'\t',line)
                if int(unit[4]) - int(unit[3]) > 100:
                    ping = '%s.%s' %(unit[0], unit[3])
                    data[ping] = 1
    return data

def gff_parse(infile, output, ping, mping_shared, mping_ref):
    ofile = open(output, 'w')
    s = re.compile(r'RIL(\d+)\_*')
    data = defaultdict(lambda : defaultdict(lambda : int()))
    with open (infile, 'r') as filehd:
        for line in filehd: 
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('#'): 
                unit  = re.split(r'\t',line)
                m     = s.search(unit[1])
                ril   = m.groups(0)[0] if m else unit[1]
                start = int(unit[3]) 
                end   = int(unit[4])
                chro  = unit[0]
                strand= unit[6]
                temp  = defaultdict(str)
                attrs = re.split(r';', unit[8])
                for attr in attrs:
                    #print attr
                    if not attr == '':
                        #print 'yes'
                        idx, value = re.split(r'\=|:', attr)
                        temp[idx] = value
                mping = '%s.%s' %(unit[0], unit[3])
                if not ping.has_key(mping):
                    data[ril]['total'] += 1
                    if mping_shared.has_key(mping):
                        data[ril]['shared'] += 1
                    elif mping_ref.has_key(mping):
                        data[ril]['ref'] += 1
                    else:
                        data[ril]['others'] += 1
                #print '%s\t%s\t%s\t%s\t%s\t%s' %(repid, chro, start, end, repname, repfam)
                #print '%s\t%s\t%s' %(chro, start, end)
    print >> ofile, 'Sample\tRef_mPing\tNB_only\tShared_NB_HEG4\tOthers'
    for ril in sorted(data.keys(), key=int):
        print >> ofile, 'RIL{}\t{}\t{}\t{}\t{}'.format(ril, data[ril]['total'], data[ril]['ref'], data[ril]['shared'], data[ril]['others'])
    ofile.close()
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--gff')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = '%s.ref_shared.table.txt' %(os.path.splitext(args.gff)[0])
 
    ping_pong    = read_gff('/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_Parental_Ping/HEG4.ALL_Filter.ping.gff')
    mping_shared = read_gff('HEG4.mping.shared.gff')
    mping_ref    = read_gff('HEG4.mping.ref_only.gff')  
    gff_parse(args.gff, args.output, ping_pong, mping_shared, mping_ref)

if __name__ == '__main__':
    main()

