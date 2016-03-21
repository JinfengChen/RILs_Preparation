#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python mPing_intergenic.py --input Somatic.mRNA.intersect

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def position(pos):
    if len(pos) > 1:
        for i in pos:
            if i == 'mRNA':
                continue
            else:
                return i
    else:
        return pos[0]

'intron length: length interval\tintron number'
def intron_len_sum(mping, total, fname):
    data = defaultdict(int) 
    inter = {
         -7:'>3501',
         -6:'3001-3500',
         -5:'3000-2500',
         -4:'2501-2000',
         -3:'2001-1500',
         -2:'1501-1000',
         -1:'1001-500',
         0:'501-0',
         1:'-1-500',
         2:'-501-1000',
         3:'-1001-1500',
         4:'-1501-2000',
         5:'-2001-2500',
         6:'-2501-3000',
         7:'-3001-3500',
         8:'-3501-4000',
         9:'-4001-4500',
         10:'-4501-5000',
         11:'-5000-5500',
         12:'<-5501'
    }
    tt = 0
    for m in sorted(mping.keys()):
        l = mping[m]
        if 1 > 0:
            tt += 1
            if l >= 3500:
                data[-7] += 1
            elif l >= 3000:
                data[-6] += 1
            elif l >= 2500:
                data[-5] += 1
            elif l >= 2000:
                data[-4] += 1
            elif l >= 1500:
                data[-3] += 1
            elif l >= 1000:
                data[-2] += 1
            elif l >= 500:
                data[-1] += 1
            elif l >= 0:
                data[0] += 1
            elif l >= -500:
                data[1] += 1
            elif l >= -1000:
                data[2] += 1
            elif l >= -1500:
                data[3] += 1
            elif l >= -2000:
                data[4] += 1
            elif l >= -2500:
                data[5] += 1
            elif l >= -3000:
                data[6] += 1
            elif l >= -3500:
                data[7] += 1
            elif l >= -4000:
                data[8] += 1
            elif l >= -4500:
                data[9] += 1
            elif l >= -5000:
                data[10] += 1
            elif l >= -5500:
                data[11] += 1
            else:
                data[12] += 1

    ofile = open(fname, 'w')
    for r in sorted(data.keys(), key=int):
        print >> ofile, '%s\t%s\t%s\t%s' %(r, inter[r], data[r], float(data[r])/int(total))


def writefile(data, fname):
    ofile = open(fname, 'w') 
    for m in sorted(data.keys()):
        print >> ofile, data[m]
    ofile.close()

'''
Chr1    RelocaTE        mPing   10046640        10046642        .       .       .       Strains=HEG4_2;GT=heterozygous  Chr1    MSU_osa1r7      mRNA    10048264        10050309        .       -       .       ID=LOC_Os01g17470.1
'''
def readtable(infile):
    data = defaultdict(list)
    dist5 = defaultdict(list)
    dist3 = defaultdict(list)
    total = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = unit[0] + '_' + unit[3] + '_' + unit[4]
                pt = unit[11]
                strand = unit[15]
                dist = int(unit[18])
                total += 1
                if dist == 0:
                    comp = 1 if abs(int(unit[3])-int(unit[12])) <= abs(int(unit[4])-int(unit[13])) else 0 # mping in gene and close to start
                    if strand == '+':
                        if comp == 1:
                            dist5[mping] = abs(int(unit[3])-int(unit[12]))
                            #print 'Five\t%s' %(line)
                        else:
                            dist3[mping] = abs(int(unit[4])-int(unit[13]))
                            #print 'Three\t%s' %(line)
                    else:
                        if comp == 1:
                            dist3[mping] = abs(int(unit[3])-int(unit[12]))
                            #print 'Three\t%s' %(line)
                        else:
                            dist5[mping] = abs(int(unit[4])-int(unit[13]))
                            #print 'Five\t%s' %(line)
                if pt == 'mRNA' and dist > 0:
                    comp = 1 if int(unit[3]) < int(unit[12]) else 0
                    if strand == '+':
                        if comp == 1:
                            dist5[mping] = -dist
                        else:
                            dist3[mping] = -dist
                    else:
                        if comp == 1:
                            dist3[mping] = -dist
                        else:
                            dist5[mping] = -dist
    prefix = os.path.splitext(infile)[0]
    intron_len_sum(dist5, total, prefix + '.5primer.distance.distr')
    intron_len_sum(dist3, total, prefix + '.3primer.distance.distr')
    writefile(dist5, prefix + '.5primer.distance.list')
    writefile(dist3, prefix + '.3primer.distance.list')
    
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

    readtable(args.input)

if __name__ == '__main__':
    main()

