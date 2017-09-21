#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def usage():
    test="name"
    message='''
Python LOD_Genome.py --input ../input/MPR.cross.uniq.QTL.mr.table

Convert position from chromosome into single genome position.

    '''
    print message

def fasta(fastafile):
    fasta_seq = defaultdict(str)
    p = re.compile(r'(\d+)')
    for record in SeqIO.parse(fastafile,"fasta"):
        m = p.search(record.id)
        rank = m.groups(0)[0] if m else 0
        fasta_seq[rank] = str(record.seq)
    return fasta_seq

'''
Convert BIN MAP
""      "GN1"   "GN10"  "GN100" "GN101" "GN102" "GN103" "GN104" "GN105" "GN106" "GN107" "GN108
"0100222046"    1       1       0       0       0       0       1       1       1       1
"0100500860"    1       1       0       0       0       0       1       1       1       1

'''
def convert_MAP(infile, outfile):
    data = defaultdict(lambda : str)
    last_pos = defaultdict(lambda : int)
    chr_end  = defaultdict(lambda : int)
    chr_end['0'] = 0
    ofile = open(outfile, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            line = line.replace('"', '')
            if line[0:1].isdigit():
                chrs = str(int(line[0:2]))
                unit = re.split(r'\t',line)
                #print '%s\t%s\t%s' %(chrs, int(unit[0][2:]), str(chr_start[chrs]))
                pos = int(unit[0][2:]) + int(chr_end[str(int(chrs)-1)])
                last_pos[chrs] = pos
                chr_end[chrs] = pos            
                #print int(unit[0][2:]), pos
                unit[0] = str(pos)
                line = '\t'.join(unit)
                print >> ofile, line
            else:
                unit = re.split(r'\t',line)
                unit[0] = 'Position'
                line = '\t'.join(unit)
                print >> ofile, line
    for c in sorted(chr_end.keys(), key=int):
        if int(c) > 0:
            c = int(c)
            mid = (int(chr_end[str(c)]) - int(chr_end[str(c-1)]))/2 + int(chr_end[str(c-1)])
            #mid = (int(chr_end[c])-int(chr_end[str(int(c)-1)]))/2 + int(chr_end[str(int(c)-1)])
            print 'Chr%s\t%s\t%s' %(c, mid, chr_end[str(c)])
    ofile.close()
    return data

def get_start(sequences):
    chr_start = defaultdict(lambda : int)
    last_start= 0
    for seq_id in sorted(sequences.keys(), key=int):
        chr_start[seq_id] = last_start
        #print seq_id, last_start, chr_start[seq_id]
        length = len(sequences[seq_id])
        last_start += length 
        midpoint = int((last_start-chr_start[seq_id])/2 + chr_start[seq_id])
        chrn = 'Chr%s' %(seq_id)
        #print '%s\t%s\t%s' %(str(seq_id), str(midpoint), str(last_start))
    return chr_start 

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

    try:
        len(args.output) > 0
    except:
        
        args.output = args.input + '.new'

    ref = '/rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa'
    fasta_seq = fasta(ref)
    #chr_start = get_start(fasta_seq)
    #print chr_start['1']
    if not os.path.isfile(args.output):
        convert_MAP(args.input, args.output)
 
if __name__ == '__main__':
    main()

