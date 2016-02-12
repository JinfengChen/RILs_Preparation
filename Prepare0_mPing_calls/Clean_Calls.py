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
python Clean_Calls.py --gff RIL275_RelocaTEi.CombinedGFF.characterized.gff
 
Clean calls with wrong TSD, or supporting junction, so every call left should have three basepair TSD, we use very strict rule to overlap.
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
                ping = '%s.%s' %(unit[0], unit[3])
                data[ping] = 1
    return data

def gff_parse(infile, output, ping):
    ofile = open(output, 'w')
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd: 
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('#'): 
                unit  = re.split(r'\t',line)
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
                if len(temp['TSD']) == 3 and not ping.has_key(mping):
                    print >> ofile, line
                #print '%s\t%s\t%s\t%s\t%s\t%s' %(repid, chro, start, end, repname, repfam)
                #print '%s\t%s\t%s' %(chro, start, end)
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
        args.output = '%s.clean.gff' %(os.path.splitext(args.gff)[0])
 
    ping = read_gff('/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_Parental_Ping/HEG4.ALL_Filter.ping.gff')
    gff_parse(args.gff, args.output, ping)

if __name__ == '__main__':
    main()

