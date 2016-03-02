#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import subprocess

def usage():
    test="name"
    message='''
python mping.excision.avg.py

Summary average number of excision events according to allele frequency in the population.
    '''
    print message

#Chr3    12409837        12409839        Chr3:12409837-12409839  +       90      0.391304347826  Parental
def readfrq(infile):
    data = defaultdict(lambda : float())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[3]] = float(unit[6])
    return data

#mPing	NumberOfEvents	NumberOfEvents_Footprint	NumberOfRILs	NumberOfRILs_Footprint	Footprint:RILs
#Chr10:11955070-11955072	1	0	1	0	Perfect:RIL185
def readfile(infile, frq):
    data = defaultdict(lambda: list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'): 
                unit = re.split(r'\t',line)
                try:
                    unit[2] = frq[unit[0]]
                except:
                    print line, frq[unit[0]]
                    unit.append(frq[unit[0]])
                unit[1] = unit[3]
                if float(unit[2]) < 0.05:
                    print line
                    data['0.025'].append(unit[1])
                elif float(unit[2]) < 0.1:
                    data['0.075'].append(unit[1])
                elif float(unit[2]) < 0.15:
                    data['0.125'].append(unit[1])
                elif float(unit[2]) < 0.2:
                    data['0.175'].append(unit[1])
                elif float(unit[2]) < 0.25:
                    data['0.225'].append(unit[1])                
                elif float(unit[2]) < 0.3:
                    data['0.275'].append(unit[1])
                elif float(unit[2]) < 0.35:
                    data['0.325'].append(unit[1])
                elif float(unit[2]) < 0.4:
                    data['0.375'].append(unit[1])
                elif float(unit[2]) < 0.45:
                    data['0.425'].append(unit[1])
                elif float(unit[2]) < 0.5:
                    data['0.475'].append(unit[1])
                elif float(unit[2]) < 0.55:
                    data['0.525'].append(unit[1])
                elif float(unit[2]) < 0.6:
                    data['0.575'].append(unit[1])
                elif float(unit[2]) < 0.65:
                    data['0.625'].append(unit[1])
                elif float(unit[2]) < 0.7:
                    data['0.675'].append(unit[1])
                else:
                    data['0.725'].append(unit[1])
    ofile = open('mping.excision.avg_frequency.table', 'w')
    for k in sorted(data.keys(), key=float):
        avg = mean(map(float ,data[str(k)]))
        var = std(map(float ,data[str(k)]))
        values = ','.join(data[str(k)])
        print >> ofile, k, avg, var, len(data[str(k)]), values
    return data
    ofile.close()

def main():
    frq = readfrq('RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.parental.frequency')
    s = readfile('Excision_newpipe_version1.footprint.list.noPing.txt', frq)

if __name__ == '__main__':
    main()




