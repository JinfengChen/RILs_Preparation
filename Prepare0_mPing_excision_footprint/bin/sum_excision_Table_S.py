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
    ofile = open('Excision_newpipe_version1.footprint.list.noPing.TableS.txt', 'w') 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'): 
                unit = re.split(r'\t',line)
                if unit[1] == '0':
                    unit[5] = 'NA'
                    unit[0] = re.split(r'-', unit[0])[0]
                    print >> ofile, '\t'.join(unit)
                    continue
                events = re.split(r';', unit[5])
                string_fp = []
                for e in sorted(events):
                    fp, rils = re.split(r':', e)
                    ril_id = re.split(r',', rils)
                    fp_short = ''
                    if '-' in fp:
                        fps = re.split(r'-', fp)
                        f_list = []
                        for f in fps:
                            f_short = '_'.join(re.split(r'_', f)[1:])
                            f_list.append(f_short)
                        fp_short = '-'.join(f_list)
                    elif '_' in fp:
                        fp_short = '_'.join(re.split(r'_', fp)[1:])
                    else:
                        fp_short = fp
                    temp = '%s(%s)' %(fp_short, str(len(ril_id)))
                    string_fp.append(temp)
                mping = re.split(r'-', unit[0])[0]
                print >> ofile, '%s\t%s\t%s\t%s\t%s\t%s' %(mping, unit[1], unit[2], unit[3], unit[4], ';'.join(string_fp))
            else:
                print >> ofile, line    
             
    ofile.close() 

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

