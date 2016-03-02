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
python Excision_Distance_Boxplot.py --highexcision --distance

--histexcision: high excision table
--distance: excision distance table from Excision_Distance.py 

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


#read gff format
def read_gff(infile, distance):
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
                mping = '%s.%s' %(chro, start)
                data[mping] = distance[mping] if distance.has_key(mping) else 'NA'
    return data

#Chr1:36267659-36267661  101,112,95,119,66,89,127,143,151
def read_highexcision(infile, distance):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'): 
                unit = re.split(r' |\t',line)
                mping = re.split(r'-', unit[0])[0]
                mping = re.sub(r':', r'.', mping)
                if distance.has_key(mping):
                    data[mping] = distance[mping]
                else:
                    data[mping] = 'NA'
                    print 'no distance : %s' %(mping)
                #data[mping] = distance[mping] if distance.has_key(mping) else 'NA'
    return data

#mPing   Excision        Distance
#Chr1.10903901   3       322234
def read_distance1(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = unit[2] 
    return data

#Chr3.29404858   Chr3.29404901   43      -       +
def read_distance(infile):
    data = defaultdict(lambda : int())
    pair = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping1_idx = unit[0]
                mping2_idx = unit[1]
                #mping1 = re.split(r'\.', unit[0])
                #mping2 = re.split(r'\.', unit[1])
                #mping1_idx = '%s_%s_%s' %(mping1[0], mping1[1], str(int(mping1[1]) + 2))
                #mping2_idx = '%s_%s_%s' %(mping2[0], mping2[1], str(int(mping2[1]) + 2))
                #mping1_idx = '%s_%s_%s' %(mping1[0], str(int(mping1[1]) + 2), mping1[1])
                #mping2_idx = '%s_%s_%s' %(mping2[0], str(int(mping2[1]) + 2), mping2[1])
                #if data.has_key(mping1_idx):
                #    print 'present in more than 1 pairs: %s' %(mping1_idx)
                #if data.has_key(mping1_idx):
                #    print 'present in more than 1 pairs: %s' %(mping2_idx)
                if data.has_key(mping1_idx):
                    data[mping1_idx] = int(unit[2]) if int(unit[2]) < data[mping1_idx] else data[mping1_idx]
                else:
                    data[mping1_idx] = int(unit[2])
                if data.has_key(mping2_idx):
                    data[mping2_idx] = int(unit[2]) if int(unit[2]) < data[mping2_idx] else data[mping2_idx]
                else:
                    data[mping2_idx] = int(unit[2])
                #data[mping1_idx] = 1
                #data[mping2_idx] = 1
                #index       = '%s:%s' %(unit[0], unit[1])
                #pair[index] = [mping1_idx, mping2_idx, unit[2], unit[3], unit[4]]
    #print 'linked mPing: %s' %(str(len(data.keys())))
    print 'mPing with distance: %s' %(len(data.keys()))
    return data


def write_distance_excision(mping, excision, distance, outfile, flag):
    ofile =open(outfile, 'w')
    print >> ofile, 'mPing\tExcision\tDistance'
    for m in sorted(mping.keys()):
        if flag == 1:
            # for excision old version we assign 0 to all mPing without value
            mping_excision = excision[m] if excision.has_key(m) else '0'
            mping_distance = distance[m] if distance.has_key(m) else 'NA'
            if int(mping_excision) < 20:
                print >> ofile, '%s\t%s\t%s' %(m, mping_excision, mping_distance)
        else:
            # for excision new version we use default
            mping_excision = excision[m] if excision.has_key(m) else 'NA'
            mping_distance = distance[m] if distance.has_key(m) else 'NA'
            if not mping_excision == 'NA':
                if int(mping_excision) < 20:
                    print >> ofile, '%s\t%s\t%s' %(m, mping_excision, mping_distance)

    ofile.close()

def writetable(highmping, allmping, prefix):
    ofile1 = open('%s.high.txt' %(prefix), 'w')
    ofile2 = open('%s.all.txt' %(prefix), 'w')
    high_dist = '\n'.join(map(str, highmping.values()))
    all_dist  = '\n'.join(map(str, allmping.values()))
    print >> ofile1,  high_dist
    print >> ofile2,  all_dist
    ofile1.close()
    ofile2.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--highexcision')
    parser.add_argument('--distance')
    parser.add_argument('--gff')
    parser.add_argument('--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.distance) > 0
    except:
        usage()
        sys.exit(2)

    prefix = 'Excision_distance.boxplot'
    if args.output:
        prefix = args.output

    distance          = read_distance(args.distance)
    highexcisionmping = read_highexcision(args.highexcision, distance)
    mping             = read_gff(args.gff, distance) 
    writetable(highexcisionmping, mping, prefix)
if __name__ == '__main__':
    main()

