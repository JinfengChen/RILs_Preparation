#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import subprocess

def usage():
    test="name"
    message='''
python footprint.py --input mping.excision.draw.example

    '''
    print message

def get_fasta_seq(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = str(record.seq)
    return fastaid

##default flanking of mping
def get_flank_seq(mping, flank_len):
    reg = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match = reg.search(mping)
    chro  = match.groups(0)[0]
    start = int(match.groups(0)[1]) - flank_len
    end   = int(match.groups(0)[2]) + flank_len
    refseq= get_fasta_seq('/rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa')
    flank = refseq[chro][start:end]
    return flank

##flanking of mping with footprint
def get_footprint_str(mping, flank, fp_dict, flank_len):
    reg = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match = reg.search(mping)
    chro  = match.groups(0)[0]
    start = int(match.groups(0)[1])
    end   = int(match.groups(0)[2])
    flank_str = list(flank)
    flank_ins = ''
    for fp_start in sorted(fp_dict.keys(), key=int):
        if fp_dict[fp_start][0] == 'D':
            for i in range(fp_start-start-(flank_len + 3), fp_start-start- (flank_len + 3) +int(fp_dict[fp_start][1])):
                flank_str[i] = '-'
        if fp_dict[fp_start][0] == 'I':
            flank_ins = ' '*(fp_start-start+ (flank_len - 3)) + '|' + fp_dict[fp_start][1]
    return [flank_ins, ''.join(flank_str)]

#Chr1:6806761-6806763	10,117,130,134,207,44,264,78,22
#Chr1:36270511-36270513	122,125,127,63,98,89
def readtable(infile, output):
    flank_len = 25 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                data = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : list()))) 
                unit = re.split(r'\t',line)
                mping = unit[0]
                flank = get_flank_seq(mping, flank_len)
                rils  = re.split(r',', unit[1])
                #print '>%s' %(mping)
                #print '%s\t%s' %('{0:15}'.format('Nipponbare'), flank)
                for ril in rils:
                    bam = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_correct_merged/GN%s.bam' %(ril)
                    subbam(mping, ril, bam, output)
                    #data.update(subbam(mping, ril, bam, output))
                    #fp_dict = data[mping][ril]
                    #flank_ril = flank
                    #fp_str  = get_footprint_str(mping, flank_ril, fp_dict, flank_len)
                    #if fp_str[0] is not '':
                    #    print '%s\t%s' %('{0:15}'.format(' '), fp_str[0])
                    #else:
                    #    print ' '
                    #print '%s\t%s' %('{0:15}'.format('RIL%s' %(ril)), fp_str[1])

def subbam(mping, ril, bam, output):
    data = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : list())))
    reg = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match = reg.search(mping)
    start = int(match.groups(0)[1])
    end   = int(match.groups(0)[2])
    chro    = match.groups(0)[0]
    region  = '%s:%s-%s' %(chro, start-100000, end+100000)
    outdir  = './%s/%s' %(output, mping)
    if not os.path.exists(output):
        os.mkdir(output)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    test_bam = '%s/%s_%s.bam' %(outdir, ril, mping)
    os.system('samtools view -hb %s %s > %s' %(bam, region, test_bam))
    os.system('samtools index %s' %(test_bam))

    cmd = 'samtools view %s %s' %(bam, mping)
    out = subprocess.check_output(cmd, shell=True)

    ##parse align around mPing
    #print 'mPing: %s' %(mping)
    #print '%s' %(cmd)
    total = 0
    footprint = []
    pattern = re.compile(r'([0-9]+)([A-Z]+)')
    lines = re.split(r'\n', out)
    for line in lines:
        unit = re.split(r'\t', line)
        if len(unit) < 2:
            continue
        #print 'Read: %s' %(unit[0])
        total += 1
        ref_start = int(unit[3])
        for (base, match) in re.findall(pattern, unit[5]):
            if match == 'M':
                ref_start += int(base)
            elif match == "I":
                ins_s   = ref_start - int(unit[3])
                ins_e   = ins_s + int(base)
                ins_seq = unit[9][ins_s:ins_e]
                footprint.append([ref_start, 'I', ins_seq])
            elif match == 'D':
                footprint.append([ref_start, 'D', base])
                ref_start += int(base)
    ##remove low frequency footprint
    temp_dict = defaultdict(lambda : int())
    fp_dict   = defaultdict(lambda : list())
    for fp in footprint:
        index = '_'.join(map(str, fp))
        temp_dict[index] += 1
    for i in range(len(footprint)):
        #print i
        index = '_'.join(map(str, footprint[i]))
        if temp_dict[index] > total*0.3:
            fp_dict[footprint[i][0]] = [footprint[i][1], footprint[i][2]]

    ##output footprint table
    for start in sorted(fp_dict.keys(), key=int): 
        #print '%s\t%s\t%s' %(start, fp_dict[start][0], fp_dict[start][1])
        data[mping][ril][start] = [fp_dict[start][0], fp_dict[start][1]] 

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

    if not args.output:
        args.output = 'test_example_draw'

    readtable(args.input, args.output)

if __name__ == '__main__':
    main()

