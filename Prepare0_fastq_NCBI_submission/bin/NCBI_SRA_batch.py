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


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Sample  #Read   Average Total   Depth   Mapped_Depth    Mapped_rate     #Library        FileName
#RIL1    23383054        100     2338305400      6.2857672043    6.13189677419   0.975520819479  1       Bam_correct/RIL1_0_CGTACG_FC153L5.recal.bam
def read_stat(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'RIL'): 
                unit = re.split(r'\t',line)
                lib  = os.path.split(unit[-1])[1]
                lib  = re.sub(r'\.recal\.bam', r'', lib)
                model = ''
                if int(unit[2]) == 100 or int(unit[2]) == 150:
                    model = 'Illumina HiSeq 2500'
                elif int(unit[2]) == 75:
                    model = 'NextSeq 500'
                elif int(unit[2]) == 250:
                    model = 'Illumina MiSeq'
                data[lib] = model
                #print lib, model
    return data

#RIL1    1       /rhome/cjinfeng/Rice/RIL/Illumina_correct/RIL1_0/RIL1_0_CGTACG_FC153L5_p1.fq.gz
def read_list(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                ril  = re.sub(r'RIL', r'', unit[0])
                files= []
                for f in unit[2:]:
                    f1 = os.path.split(f)[1]
                    #f2 = re.sub(r'p1', r'p2', f1)
                    files.append(f1)
                    #files.append(f2)
                data[ril] = files
    return data

def NCBI_SRA(rils, model, tsv):
    ofile = open(tsv, 'a')
    for ril in sorted(rils.keys(), key=int):
        #print ril
        #print >> ofile, 'RIL%s\tnot collected\tPRJNA316308\tOryza Sativa\tnot collected\tNipponbare x HEG4 RIL%s\tTemperate Japonica\tnot collected\tSeedling\tnot collected\tLeaf\tSusan Wessler\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected\tnot collected' %(ril, ril)
        #RIL1_0_CGTACG_FC153L5_p1.fq.gz
        for lib in rils[ril]:
            lib_id   = re.sub(r'_p1.fq.gz', r'', lib)
            pl_model = 'Illumina HiSeq 2500' if not model.has_key(lib_id) else model[lib_id]
            lib_files= '%s\t%s' %(lib, re.sub(r'p1', r'p2', lib))
            print >> ofile, 'RIL%s\tPRJNA316308\tRIL%s\t%s\tnot collected\tWGS\tGENOMIC\tRANDOM\tPaired\tILLUMINA\t%s\tfastq\t%s' %(ril, ril, lib_id, pl_model, lib_files)
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

    #cp SRA_metadata.tsv RILs_272.NCBI_SRA.tsv
    os.system('cp ../input/SRA_metadata.tsv RILs_272.NCBI_SRA.tsv')
    rils  = read_list(args.input)
    model = read_stat('Bam_correct.bam.stat') 
    NCBI_SRA(rils, model, 'RILs_272.NCBI_SRA.tsv')

if __name__ == '__main__':
    main()

