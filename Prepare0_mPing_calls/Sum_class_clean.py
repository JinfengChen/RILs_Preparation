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
python Sum_class_clean.py --input RIL275_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff

classify mPing into parental, shared in ril, or unique and then classify into hom, het, som
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def get_class(infile):
    data = defaultdict(lambda : int())
    r1 = re.compile(r'type=hom')
    r2 = re.compile(r'type=het')
    r3 = re.compile(r'type=som')
    with open (infile, 'r') as filehd:
        for line in filehd: 
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('#'):
                if r1.search(line):
                    data[0] += 1
                    data[3] += 1
                elif r2.search(line):
                    data[1] += 1
                    data[3] += 1
                elif r3.search(line):
                    data[2] += 1
                    data[3] += 1
    return data

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
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    prefix = os.path.splitext(os.path.splitext(args.input)[0])[0]
    os.system('bedtools intersect -a %s.gff -b HEG4.ALL.mping.non-ref.AF0.1.gff > %s.shared_parental.gff' %(prefix, prefix))
    os.system('bedtools intersect -a %s.gff -b HEG4.ALL.mping.non-ref.AF0.1.gff -v > %s.nonparental.gff' %(prefix, prefix))
    os.system('bedtools intersect -a %s.nonparental.gff -b %s.unique_mPing.gff -v > %s.shared_nonparental.gff' %(prefix, prefix, prefix))
 
    unique_class = get_class(args.input)
    shared_non_class = get_class('%s.shared_nonparental.gff' %(prefix)) 
    shared_p_class   = get_class('%s.shared_parental.gff' %(prefix))
    
    ofile  = open ('%s.class.summary' %(prefix), 'w') 
    print >> ofile, 'Type\tHomozygous\tHeterzygous\tSomatic'
    #print >> ofile, 'Parental\t%s\t%s\t%s' %(float(shared_p_class[0])/shared_p_class[3], float(shared_p_class[1])/shared_p_class[3], float(shared_p_class[2])/shared_p_class[3])
    #print >> ofile, 'Unique\t%s\t%s\t%s' %(float(unique_class[0])/unique_class[3], float(unique_class[1])/unique_class[3], float(unique_class[2])/unique_class[3])
    #print >> ofile, 'Shared\t%s\t%s\t%s' %(float(shared_non_class[0])/shared_non_class[3], float(shared_non_class[1])/shared_non_class[3], float(shared_non_class[2])/shared_non_class[3])
    print >> ofile, 'Parental\t%s\t%s\t%s' %(float(shared_p_class[0]), float(shared_p_class[1]), float(shared_p_class[2]))
    print >> ofile, 'Shared\t%s\t%s\t%s' %(float(shared_non_class[0]), float(shared_non_class[1]), float(shared_non_class[2]))
    print >> ofile, 'Unique\t%s\t%s\t%s' %(float(unique_class[0]), float(unique_class[1]), float(unique_class[2]))
    ofile.close()

if __name__ == '__main__':
    main()

