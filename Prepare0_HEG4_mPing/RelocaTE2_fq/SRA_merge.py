#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob
import gzip

def usage():
    test="name"
    message='''
python ReNameSRA_merge.py --input ./CAAS_1/ --ouput CAAS_1_fastq

Read all directory under CAAS_1, convert sra into fastq and merge fastq in one subdir into two _1 and _2 files
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 30 --lines %s --interval 120 --task 12 --mem 20G --time 10:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

#ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP003/SRP003189/
#"SRX025294","Resequencing of 50 rice individuals-rufipogon_Yuan3-9 ","Oryza sativa","Illumina Genome Analyzer II","BGI","SRP003189","Resequencing of 50 rice individuals","SRS086373","","243.58","1","6711790","590637520","ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX025/SRX025294","SZC08001CTDCAAPE","WGS","GENOMIC","PCR"
def read_sra_inf(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit    = re.split(r',', line)
                acc     = re.sub(r'"','', unit[0])
                strain  = re.sub(r'"','', unit[1][37:])
                strain  = re.sub(r' ','', strain)
                link    = re.sub(r'"','', unit[13])
                link    = re.sub(r' ','', link)
                data[acc] = [strain, link]
    return data


#SRR063622,2011-12-15 10:40:33,2011-12-15 10:40:33,33451859,6690371800,33451859,200,3134,,ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR063/SRR063622/SRR063622.sra,SRX025244
def read_sra_run(infile, sra2strain):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit    = re.split(r',', line)
                acc     = unit[0]
                run     = unit[10]
                if sra2strain.has_key(run):
                    #print acc, sra2strain[run]
                    data[acc] = sra2strain[run]
    return data


def merge_fq(fastq, fq_out):
    infile = gzip.open(fastq, 'r')
    ofile  = gzip.open(fq_out, 'ab')
    for record in SeqIO.parse(infile, "fastq"):
        SeqIO.write(record, ofile, 'fastq')
    infile.close()
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
    
    
    strain2fq = defaultdict(lambda : defaultdict(lambda : list()))
    sras = glob.glob('%s/*/*.sra' %(args.input))
    ofile = open('dump.sh', 'w')
    for sra in sras:
        dirname= os.path.abspath(os.path.dirname(sra))
        acc    = os.path.split(os.path.splitext(sra)[0])[-1] 
        cmd    = []
        fq1    = '%s/%s_1.fastq' %(dirname, acc)
        fq2    = '%s/%s_2.fastq' %(dirname, acc)
        fq1gz  = '%s/%s_1.fastq.gz' %(dirname, acc)
        fq2gz  = '%s/%s_2.fastq.gz' %(dirname, acc)
        if not os.path.isfile(fq1) and not os.path.isfile(fq1gz):
            cmd.append('/opt/linux/centos/7.x/x86_64/pkgs/sratoolkit/2.5.0/bin/fastq-dump --outdir %s --split-files %s' %(dirname, os.path.abspath(sra)))
            cmd.append('pigz -p 12 %s' %(fq1))
            cmd.append('pigz -p 12 %s' %(fq2))
            #pass
        elif os.path.isfile(fq1) and not os.path.isfile(fq1gz):
            cmd.append('pigz -p 12 %s' %(fq1))
            cmd.append('pigz -p 12 %s' %(fq2))
        if len(cmd) > 1:
            print >> ofile, '\n'.join(cmd)
    ofile.close()
    runjob('dump.sh', 3) 

if __name__ == '__main__':
    main()

