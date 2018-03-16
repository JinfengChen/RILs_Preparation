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


def runjob(script, lines, cpu, queue):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines %s --interval 120 --task %s --mem 15G --time 10:00:00 --queue %s --convert no %s' %(lines, cpu, queue, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr5:22293155-22293157  RIL     270,222
def read_list(infile):
    data = defaultdict(lambda : defaultdict(lambda : str()))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                rils = re.split(r',', unit[2])
                for ril in rils:
                    data[ril][unit[0]] = unit[1]
    return data

def sort_mping(mping_list, infile):
    list_up    = []
    list_down  = []
    list_all   = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                if not mping_list.has_key(unit[0]):
                    continue
                rils = re.split(r',', unit[2])
                up   = 0
                down = 0
                for ril in rils:
                    if int(ril) <= 139:
                        up += 1
                    else:
                        down += 1
                if up > 0 and down == 0:
                    list_up.append(unit[0]) 
                elif up > 0 and down > 0 and float(up)/float(down) > 2.0:
                    list_up.append(unit[0])
                elif up == 0 and down > 0:
                    list_down.append(unit[0])
                elif up > 0 and down > 0 and float(up)/float(down) < 0.5:
                    list_down.append(unit[0])
                else:
                    list_all.append(unit[0])
                print '%s\t%s\t%s' %(unit[0], up, down)
    list_final = []
    list_final.extend(list_up)
    list_final.extend(list_down)
    list_final.extend(list_all)
    print '\t'.join(list_up)
    print '\t'.join(list_down)
    print '\t'.join(list_all)
    print '\t'.join(list_final)
    return list_final

def readtable(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 0:
                unit = re.split(r'\t',line)
                data[unit[0]] = 1
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

    parental_list = readtable('Shared_Parental_0.3.list')
    shared_list   = readtable('Shared_RIL_0.1.list')
    rils          = readtable('RILs.list')
    mping_mat     = read_list(args.input)
    parental_list_sorted = sort_mping(parental_list, args.input)
    shared_list_sorted   = sort_mping(shared_list, args.input)

    ofile_p = open('Shared_Parental_0.3.matrix', 'w')
    ofile_s = open('Shared_RIL_0.1.matrix', 'w')
    print >> ofile_p, '%s' %('\t'.join(parental_list_sorted))
    print >> ofile_s, '%s' %('\t'.join(shared_list_sorted))
    for ril in sorted(rils.keys(), key=int):
        score = []
        #temp  = defaultdict(lambda : str())
        for mping_p in parental_list_sorted:
            flag_p = '0'
            if mping_mat.has_key(ril):
                if mping_mat[ril].has_key(mping_p):
                    flag_p = '1'
            #temp[mping_p] = flag_p
            score.append(flag_p)
        print >> ofile_p, '%s' %('\t'.join(score))
        score = []
        for mping_s in shared_list_sorted:
            flag_s = '0'
            if mping_mat.has_key(ril):
                if mping_mat[ril].has_key(mping_s):
                    flag_s = '1'
            score.append(flag_s)
        print >> ofile_s, '%s' %('\t'.join(score))
    ofile_p.close()
    ofile_s.close()

if __name__ == '__main__':
    main()

