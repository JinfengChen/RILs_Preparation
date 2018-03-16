#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python Sum_unique.py --input RIL275_RelocaTEi.CombinedGFF.characterized.unique_mPing.gff

get list of unique mPing number for each RILs, generate mping frquency for shared mPing in RILs, can check parental mPing, shared in RILs
 
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#split unit[8] of gff
def gff_attr8(unit8):
    temp = defaultdict(lambda : str())
    attrs = re.split(r';', unit8)
    for attr in attrs:
        #print attr
        if not attr == '':
            #print 'yes'
            idx, value = re.split(r'\=', attr)
            temp[idx] = value
    return temp

#correct mping index for these not accurate calls
def readtable_ril_mping_correct(infile):
    mping_correct = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not unit[1] == unit[10]:
                    mping= '%s:%s-%s' %(unit[0], unit[3], unit[4])
                    attrs1 = gff_attr8(unit[8])
                    attrs2 = gff_attr8(unit[17])
	            if attrs1['TSD'] == 'supporting_junction' or not len(attrs1['TSD']) == 3:
                        if not mping_correct.has_key(mping):
                            if not attrs2['TSD'] == 'supporting_junction' and len(attrs2['TSD']) == 3:
                                mping_correct[mping] = '%s:%s-%s' %(unit[9], unit[12], unit[13])
    return mping_correct


##overlap with ril
#Chr1    RIL231_0        transposable_element_attribute  4228091 4228092 +       .       .       ID=Chr1.4228092.spanners;Strain=RIL231_0;avg_flankers=6;spanners=0;type=homozygous;TE=mping;TSD=TT      Chr1    RIL231_0
#some of mPing insertion sites are not accurate. we create a dict to store correct index of this mping using their overlap.
#the resulted allele frequency should have correct position for all the mping
def readtable_ril(infile, mping_correct):
    data = defaultdict(lambda : defaultdict(lambda : int()))
    r = re.compile(r'RIL(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                #if not unit[1] == unit[10]:
                if not unit[1] == unit[10] and unit[3] == unit[12] and unit[4] == unit[13]:
                    mping1= '%s:%s-%s' %(unit[0], unit[3], unit[4])
                    mping2= '%s:%s-%s' %(unit[9], unit[12], unit[13])
                    ril1  = r.search(unit[1]).groups(0)[0] if r.search(unit[1]) else 'NA'
                    ril2  = r.search(unit[10]).groups(0)[0] if r.search(unit[10]) else 'NA'
                    if mping_correct.has_key(mping1):
                        mping1 = mping_correct[mping1]
                    if mping_correct.has_key(mping2):
                        mping2 = mping_correct[mping2]
                    #print '%s\t%s\t%s\t%s' %(mping1, ril1, mping2, ril2)
                    data[mping1][ril1] = 1
                    data[mping2][ril2] = 1
    return data

##overlap with HEG4
def readtable_nonref(infile, mping_correct):
    data = defaultdict(lambda : defaultdict(lambda : int()))
    r = re.compile(r'RIL(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                if unit[3] == unit[12] and unit[4] == unit[13]:
                    mping= '%s:%s-%s' %(unit[0], unit[3], unit[4])
                    ril  = r.search(unit[1]).groups(0)[0] if r.search(unit[1]) else 'NA'
                    if mping_correct.has_key(mping):
                        mping = mping_correct[mping]
                    data[mping][ril] = 1
    return data


##unique mping
def readtable(infile):
    data = defaultdict(lambda : int())
    r = re.compile(r'RIL(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                ril  = r.search(unit[1]).groups(0)[0] if r.search(unit[1]) else 'NA'
                data[ril] += 1
    return data

##get class for each ril from gff
#Chr1    RIL231_0        transposable_element_attribute  36023407        36023409
def get_ril_class(infile):
    data = defaultdict(lambda : defaultdict(lambda : int()))
    r1 = re.compile(r'type=hom')
    r2 = re.compile(r'type=het')
    r3 = re.compile(r'type=som')
    with open (infile, 'r') as filehd:
        for line in filehd: 
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('#'):
                unit = re.split(r'\t',line)
                ril  = re.split(r'_',unit[1])[0]
                ril  = re.sub(r'RIL', r'', ril)
                if r1.search(line):
                    data[ril][0] += 1
                    data[ril][3] += 1
                elif r2.search(line):
                    data[ril][1] += 1
                    data[ril][3] += 1
                elif r3.search(line):
                    data[ril][2] += 1
                    data[ril][3] += 1
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

    prefix = os.path.splitext(os.path.splitext(args.input)[0])[0]

    #mping_correct_index  = defaultdict(lambda : str()) 
    mping_correct_index = readtable_ril_mping_correct('%s.overlap_ril' %(prefix))
    mping_ovlp_rils          = readtable_ril('%s.overlap_ril' %(prefix), mping_correct_index)
    mping_ovlp_heg4          = readtable_nonref('%s.overlap_ref' %(prefix), mping_correct_index)
    mping_ovlp_nb_rils           = readtable_nonref('%s.overlap_ref_ril' %(prefix), mping_correct_index)
    mping_ovlp_nb_ref            = readtable_nonref('%s.overlap_ref_NB' %(prefix), mping_correct_index)

    r = re.compile(r'(\w+):(\d+)-(\d+)')
    ##mPing_allele_frequency
    ofile = open('%s.shared_mping.ril.frequency' %(prefix), 'w') 
    ofile1 = open('%s.shared_mping.ril.list' %(prefix), 'w')
    for mping in mping_ovlp_rils.keys():
        m = r.search(mping)
        chro, start, end = ['', 0, 0]
        if m:
            chro = m.groups(0)[0]
            start = m.groups(0)[1]
            end   = m.groups(0)[2]
        count = len(mping_ovlp_rils[mping].keys())
        if mping_ovlp_heg4.has_key(mping):
            print >> ofile1, '%s\tParental\t%s' %(mping, ','.join(map(str, mping_ovlp_rils[mping].keys())))
            print >> ofile, '%s\t%s\t%s\t%s\t%s\t%s\t%s\tParental' %(chro, start, end, mping, '+', count, float(count)/272)
        else:
            print >> ofile1, '%s\tRIL\t%s' %(mping, ','.join(map(str, mping_ovlp_rils[mping].keys())))
            print >> ofile, '%s\t%s\t%s\t%s\t%s\t%s\t%s\tRIL' %(chro, start, end, mping, '+', count, float(count)/272)
    for mping in mping_ovlp_nb_rils.keys():
        m = r.search(mping)
        chro, start, end = ['', 0, 0]
        if m:
            chro = m.groups(0)[0]
            start = m.groups(0)[1]
            end   = m.groups(0)[2]
        count = len(mping_ovlp_nb_rils[mping].keys())
        if mping_ovlp_nb_ref.has_key(mping): 
            print >> ofile1, '%s\tParental\t%s' %(mping, ','.join(map(str, mping_ovlp_nb_rils[mping].keys())))
            print >> ofile, '%s\t%s\t%s\t%s\t%s\t%s\t%s\tParental' %(chro, start, end, mping, '+', count, float(count)/272)
        else:
            print >> ofile1, '%s\tRIL\t%s' %(mping, ','.join(map(str, mping_ovlp_nb_rils[mping].keys())))
            print >> ofile, '%s\t%s\t%s\t%s\t%s\t%s\t%s\tRIL' %(chro, start, end, mping, '+', count, float(count)/272)
        
    ofile.close()
    ofile1.close()

    ##RILs shared and unique mPing
    #shared with ril
    ril_mping_count = defaultdict(lambda : int())
    for mping in mping_ovlp_rils.keys():
        for ril in mping_ovlp_rils[mping].keys():
            if mping_ovlp_heg4[mping][ril] == 0:
                ril_mping_count[ril] += 1
    #shared with heg4
    heg4_mping_count = defaultdict(lambda : int())
    for mping in mping_ovlp_heg4.keys():
        for ril in mping_ovlp_heg4[mping].keys():
            if mping_ovlp_heg4[mping][ril] > 0:
                heg4_mping_count[ril] += 1
 
    #unique 
    unique_mping = readtable(args.input)
    unique_mping_class = get_ril_class(args.input)
    #output table
    ofile = open('%s.mping.shared_unique_table.txt' %(prefix), 'w')
    print >> ofile, 'Sample\tShared_HEG4\tShared_RILs\tShared\tUnique\tUnique_hom\tUnique_het\tUnique_som'
    for ril in sorted(heg4_mping_count.keys(), key=int):
        shared_heg4 = heg4_mping_count[ril]
        shared_rils = ril_mping_count[ril]
        shared      = int(shared_heg4) + int(shared_rils)
        unique      = unique_mping[ril]
        unique_hom  = unique_mping_class[ril][0]
        unique_het  = unique_mping_class[ril][1]
        unique_som  = unique_mping_class[ril][2]
        print >> ofile, 'RIL%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(ril, shared_heg4, shared_rils, shared, unique, unique_hom, unique_het, unique_som)
    ofile.close()

    print 'Sample\tUnique_mPing'
    #unique_mping = readtable(args.input)
    for ril in sorted(unique_mping.keys(), key=int):
        print 'RIL%s\t%s' %(ril, unique_mping[ril])

if __name__ == '__main__':
    main()

