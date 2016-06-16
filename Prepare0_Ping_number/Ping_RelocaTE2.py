#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
import glob
import time
from Bio import SeqIO
import fnmatch

def usage():
    test="name"
    message='''
python ReNameSRA_sum_Ping.py --input Other_fastq_RelocaTEi_Ping

Summary RelocatTE call in current direcory using simulation in input directory
--call: RelocaTE, TEMP or other
--check: check unfinished job without results directory or non_ref.gff
--input: dir of simulation, where we can find insertion simulated gff file "MSU7.Chr4.mPing.rep1.gff"

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


#ping	ping	ping_code
#Chr1:6806761-6806763	HEG4	2
#Chr1:36270511-36270513	NB	3
def read_mping_list(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = re.split(r'-', unit[0])[0]
                data[mping] = unit[3]
                #print mping, rils
    return data

def parse_gff(infile):
    num = 0
    if not os.path.isfile(infile) or not os.path.getsize(infile):
        return 0 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                num += 1 
    return num

#Chr3    not.give        RelocaTE_i      28019800        28019802        .       +       .       ID=repeat_Chr3_28019800_28019802;
def parse_gff_ping(infile, ping_inf):
    num = 0
    if not os.path.isfile(infile) or not os.path.getsize(infile):
        return 0, [''], ['']
    ping_number = 0
    ping_code   = []
    ping_loci   = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                unit = re.split(r'\t', line)
                ping_id = '%s:%s' %(unit[0], unit[3])
                if ping_inf.has_key(ping_id):
                    ping_code.append(ping_inf[ping_id])
                else:
                    ping_code.append('*')
                ping_number += 1
                ping_loci.append(ping_id)
    return ping_number, ping_code, ping_loci

#1       B001    China   Temperate japonica      ERS470219       anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/ERR/ERR622/ERR622583/ERR622583.sra
def parse_inf(infile, sufix):
    data = defaultdict(lambda : str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                unit = re.split(r'\t', line)
                index = '%s_%s' %(unit[4], sufix)
                if not data.has_key(index):
                    data[index]= '%s\t%s\t%s' %(unit[1], unit[2], unit[3])
    return data



def parse_overlap_reloacte(infile, ref_te, call_te):
    #Total number of call, true call, call with breakpoint near TSD, No call, False call
    data = defaultdict(lambda : int)
    dupli= defaultdict(lambda : int())
    true = 0
    tsd  = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                #print line
                if dupli.has_key('%s_%s' %(unit[12], unit[13])):
                    continue
                else:
                    dupli['%s_%s' %(unit[12], unit[13])] == 1
                #print 'pass'
                pos  = map(int, [unit[3], unit[4], unit[12], unit[13]])
                dist_min = min([abs(pos[2]-pos[0]), abs(pos[2]-pos[1]), abs(pos[3]-pos[0]), abs(pos[3]-pos[0])])
                dist_max = max([abs(pos[2]-pos[0]), abs(pos[2]-pos[1]), abs(pos[3]-pos[0]), abs(pos[3]-pos[0])])
                true += 1
                if dist_min <= 5 and dist_max <= 5:
                    tsd += 1
    call = len(call_te.keys())
    ref  = len(ref_te.keys())
    data = [call, true, tsd, ref-true, call-true]
    return data




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input') 
    parser.add_argument('-l', '--list')
    parser.add_argument('-c', '--call')
    parser.add_argument('-ck', '--check')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.check:
        args.check = 0

    if not args.call: 
        args.call = 'RelocaTEi'


    ping_inf = read_mping_list('Ping.list')

    #inf  = parse_inf(args.list, args.call)
    inf = defaultdict(lambda : str())
    data = defaultdict(lambda : defaultdict(lambda : list)) 
    project = os.path.split(args.input)[1]
    sum_file = '%s.summary' %(args.input)
    ofile = open(sum_file, 'w')
    #print >> ofile, 'Accession\tPing\tPing_NonRef\tPing_Ref\tName\tOrigin\tClass'
    print >> ofile, 'Accession\tPing\tPing_code\tPing_loci'
    for call in sorted(os.listdir(args.input)):
        if not fnmatch.fnmatch(call, '*_RelocaTEi'):
            continue
        dirname = os.path.abspath('%s/%s' %(args.input, call))
        #print dirname
        if args.call == 'RelocaTE':
            pass
        elif args.call == 'RelocaTEi':
            #cmd_ping    = 'python PickPing.py --input %s --SNP' %(dirname)
            ping_nonref_gff  = '%s/repeat/results/ALL.all_nonref_insert.Ping.gff' %(dirname)
            ping_ref_gff     = '%s/repeat/results/ALL.all_ref_insert.Ping.gff' %(dirname)
            #if not os.path.exists(ping_gff):
            #os.system(cmd_ping) 
            #non_ref_gff = '%s/repeat/results/ALL.all_nonref_insert.gff' %(dirname)
            #ref_gff     = '%s/repeat/results/ALL.all_ref_insert.gff' %(dirname)
            call_nonref_ping_number, call_nonref_ping_code, call_nonref_ping_loci  = parse_gff_ping(ping_nonref_gff, ping_inf)
            call_ref_ping_number, call_ref_ping_code, call_ref_ping_loci           = parse_gff_ping(ping_ref_gff, ping_inf)
            call_ping_number  = call_nonref_ping_number + call_ref_ping_number
            call_ping_code    = call_nonref_ping_code + call_ref_ping_code
            call_ping_loci    = call_nonref_ping_loci + call_ref_ping_loci
            if call_ping_number > 0:
                print >> ofile, '%s\t%s\t%s\t%s' %(call, call_ping_number, ''.join(sorted(map(str, call_ping_code))), ';'.join(sorted(call_ping_loci)))
            else:
                print >> ofile, '%s\t0\tNA\tNA' %(call)
            #call_nonref = parse_gff(non_ref_gff)
            #call_ref    = parse_gff(ref_gff)
            #call_inf    = inf[call] if inf.has_key(call) else 'NA\tNA\tNA' 
            #if int(args.check) == 1 and (not os.path.exists('%s/repeat/results' %(dirname)) or not os.path.getsize(non_ref_gff) > 0):
            #    print >> ofile, '%s\t%s\t%s\t%s\t%s\tNo output' %(call, call_ping, call_ref, call_nonref, call_inf)
            #print >> ofile, '%s\t%s\t%s\t%s\t%s' %(call, call_ping, call_nonref_ping, call_ref_ping, call_inf)
        elif args.call == 'TEMP':
            pass
    ofile.close()  

if __name__ == '__main__':
    main()

