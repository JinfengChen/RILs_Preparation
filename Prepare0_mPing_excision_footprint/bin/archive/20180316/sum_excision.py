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
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#mPing   NumberOfEvents  NumberOfEvents_Footprint        NumberOfRILs    NumberOfRILs_Footprint  Footprint:RILs
#Chr10:11955070-11955072 1       0       1       0       Perfect:RIL185
#Chr10:16851282-16851284 3       2       4       2       16851260_D_20:RIL205;16851284_D_1-16851292_D_1-16851297_D_3-16851304_D_1:RIL237;Perfect:RIL92,RIL270
def readtable(infile):
    data = defaultdict(str)
    num_event = 0
    num_ril   = 0
    num_event_fp = 0
    num_ril_fp = 0
    num_mping_excision = 0
    num_mping_fp = 0
    ofile = open('Excision_newpipe_version1.footprint.mPing.RILs_footprint.list.txt', 'w')
   
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'): 
                unit = re.split(r'\t',line)
                if int(unit[1]) == 0:
                    continue
                print >> ofile, '>%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], unit[3], unit[4])
                num_event += int(unit[1])
                num_event_fp += int(unit[2])
                num_ril   += int(unit[3])
                num_ril_fp += int(unit[4])
                num_mping_excision += 1 if int(unit[1]) > 0 else 0
                num_mping_fp += 1 if int(unit[2]) > 0 else 0
                events = re.split(r';', unit[5])
                dict_t = defaultdict(lambda : str())
                for e in sorted(events):
                    fp, rils = re.split(r':', e)
                    rils = re.sub(r',', r'/', rils)
                    dict_t[rils] = fp
                for r in sorted(dict_t.keys()):
                    print >> ofile, '%s\t%s' %(r, dict_t[r])

    ofile.close() 
    print 'Number of excisions: %s' %(num_ril)
    print 'Number of excisions with footprint: %s;%s' %(num_ril_fp, float(num_ril_fp)/num_ril)
    print 'Number of events: %s' %(num_event)
    print 'Number of events with footprints: %s;%s' %(num_event_fp, float(num_event_fp)/num_event)
    print 'Number of mPing excised: %s' %(num_mping_excision)
    print 'Number of mPing leave footprint: %s;%s' %(num_mping_fp, float(num_mping_fp)/num_mping_excision)

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

