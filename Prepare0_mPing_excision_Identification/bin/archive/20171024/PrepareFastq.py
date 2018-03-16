#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob

def usage():
    test="name"
    message='''
python PrepareFastq.py --bam RIL_ALL_bam --gff 

Get unmapped fastq and fastq that within 100 kb of mPing. Because some mPing are linked, these will cause some reads are present more than once. We need to deal with the mPing position first, merge 10 kb interval of mPing. Use that to extract reads from bam file.

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=1,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)


#Chr1    1032974 1232977
def readtable(infile):
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                idx  = '%s:%s-%s' %(unit[0], unit[1], unit[2])
                data.append(idx)
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.bam) > 0
    except:
        usage()
        sys.exit(2)

    bam2fastq='/rhome/cjinfeng/BigData/software/bam2fastq/bam2fastq-1.1.0/bam2fastq'
    samtools ='/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools'
    if not os.path.exists('../input/Parent.ALL.mPing.100kb_flank.merge.table'):
        os.system('bedtools slop -i ../input/Parent.ALL.mPing.gff -g /rhome/cjinfeng/BigData/00.RD/seqlib/MSU7.chrlen -b 100000 > ../input/Parent.ALL.mPing.100kb_flank.gff')
        os.system('bedtools merge -i ../input/Parent.ALL.mPing.100kb_flank.gff > ../input/Parent.ALL.mPing.100kb_flank.merge.table')
    bam0s = glob.glob('%s/*.bam' %(args.bam))

    #output directory
    outdir = os.path.abspath('../input/RILs_ALL_unmapped_mping_fastq')
    if not os.path.exists(outdir):
        os.mkdir(outdir) 
    #mping region
    mping_regs = os.path.abspath('../input/Parent.ALL.mPing.100kb_flank.merge.table')

    cmd = []
    for bam in sorted(bam0s):
        bam      = os.path.abspath(bam)
        ril      = os.path.split(bam)[1]
        ril      = re.sub(r'.bam', r'', ril)
        ril      = re.sub(r'GN', r'RIL', ril)
        ril_dir  = os.path.abspath(ril)
        if not os.path.exists(ril_dir):
            os.mkdir(ril_dir)
        #unmapped bam
        #An unmapped read whose mate is mapped. 4 itself is unmapped, 264 mate is unmapped
        cmd.append('%s view -hb -f 4 -F264 %s > %s/%s.unmapped1.bam' %(samtools, bam, ril_dir, ril))
        #A mapped read who's mate is unmapped. 8 mate is unmapped, 260 itself is unmapped
        cmd.append('%s view -hb -f 8 -F260 %s > %s/%s.unmapped2.bam' %(samtools, bam, ril_dir, ril))
        #Both reads of the pair are unmapped. 12 reads are paired, 256 not primary alignment
        cmd.append('%s view -hb -f 12 -F256  %s > %s/%s.unmapped3.bam' %(samtools, bam, ril_dir, ril))
        #merge
        cmd.append('%s merge %s/%s.unmapped.bam %s/%s.unmapped[123].bam' %(samtools, ril_dir, ril, ril_dir, ril))
        #bam2fastq
        cmd.append('%s %s/%s.unmapped.bam -o %s/%s.unmapped#.fq' %(bam2fastq, ril_dir, ril, ril_dir, ril))
        
        #mping regions
        cmd.append('%s view -hb -L %s %s > %s/%s.mping.bam' %(samtools, mping_regs, bam, ril_dir, ril))
        #filter unmapped
        cmd.append('%s view -hb -f 3 %s/%s.mping.bam > %s/%s.mping.mapped.bam' %(samtools, ril_dir, ril, ril_dir, ril))
        #bam2fastq
        cmd.append('%s %s/%s.mping.mapped.bam -o %s/%s.mping.mapped#.fq' %(bam2fastq, ril_dir, ril, ril_dir, ril))   
         
        #merge fastq and mv
        cmd.append('cat %s/%s.unmapped_1.fq %s/%s.mping.mapped_1.fq > %s/%s_1.fq' %(ril_dir, ril, ril_dir, ril, ril_dir, ril))
        cmd.append('cat %s/%s.unmapped_2.fq %s/%s.mping.mapped_2.fq > %s/%s_2.fq' %(ril_dir, ril, ril_dir, ril, ril_dir, ril))
        cmd.append('rm %s/%s.unmapped* %s/%s.mping*' %(ril_dir, ril, ril_dir, ril))
        cmd.append('mv %s %s' %(ril_dir, outdir))

    ofile = open('prepare_reads.sh', 'w')
    print >> ofile, '\n'.join(cmd)
    ofile.close()

    runjob('prepare_reads.sh', 12)
if __name__ == '__main__':
    main()

