#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import subprocess
sys.path.append('%s/lib' %(os.getcwd()))

def usage():
    test="name"
    message='''
python Ping_genotype.py --input Ping.list
python Ping_genotype.py --input Ping.list --snp_map RILs_275/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam --output RIL275_ping_genotype

input files:
--input Ping.list
--snp_map: snp map to genotype ping
--bam: bam files where we can get the link information to extract library name

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def binarySearch(data, val):
    highIndex = len(data)-1
    lowIndex = 0
    while highIndex > lowIndex:
            index = (highIndex + lowIndex) / 2
            sub = int(data[index])
            #print highIndex, index, lowIndex, sub, val
            if data[lowIndex] == val:
                    return [lowIndex, lowIndex]
            elif sub == val:
                    return [index, index]
            elif data[highIndex] == val:
                    return [highIndex, highIndex]
            elif sub > val:
                    if highIndex == index:
                            return sorted([highIndex, lowIndex])
                    highIndex = index
            else:
                    if lowIndex == index:
                            return sorted([highIndex, lowIndex])
                    lowIndex = index
    return sorted([highIndex, lowIndex])

#Convert BIN MAP
'''
""      "GN1"   "GN10"  "GN100" "GN101" "GN102" "GN103" "GN104" "GN105" "GN106" "GN107" "GN108
"0100222046"    1       1       0       0       0       0       1       1       1       1
"0100500860"    1       1       0       0       0       0       1       1       1       1
'''

def convert_MAP(infile):
    rils = []
    data = defaultdict(lambda : defaultdict(lambda: defaultdict(lambda: str)))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            line = line.replace('"', '')
            if line[0:1].isdigit():
                unit = re.split(r'\t',line)
                #print '%s\t%s\t%s' %(chrs, int(unit[0][2:]), str(chr_start[chrs]))
                chrs = 'Chr%s' %(str(int(unit[0][0:2])))
                pos = int(unit[0][2:])
                for i in range(1,len(unit)):
                    #print i, rils[i], chrs, pos
                    ril = rils[i]
                    data[ril][chrs][pos] = unit[i]
            else:
                unit = re.split(r'\t',line)
                unit[0] = 'Position'
                #print unit[0], unit[1], unit[2]
                rils.extend(unit)
   
    #for t in sorted(data['GN204']['Chr10'].keys(), key=int):
    #    print t
    return data

#Convert SNP map
#""      "GN1"   "GN10"  "GN100" "GN101" "GN102" "GN103" "GN104" "GN105" "GN106" 
#"0100021547A"   NA      NA      NA      0       0       0       NA      NA      
#"0100031071A"   NA      1       0       0       0       0       1       1       
#"0100031478C"   1       1       0       0       0       0       NA      1       

def convert_MAP_SNP(infile):
    rils = []
    data = defaultdict(lambda : defaultdict(lambda: defaultdict(lambda: str)))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            line = line.replace('"', '')
            if line[0:1].isdigit():
                unit = re.split(r'\t',line)
                #print '%s\t%s\t%s' %(chrs, int(unit[0][2:]), str(chr_start[chrs]))
                chrs = 'Chr%s' %(str(int(unit[0][0:2])))
                unit[0] = unit[0][:-1] #remove reference base
                pos = int(unit[0][2:])
                for i in range(1,len(unit)):
                    #print i, rils[i], chrs, pos, unit[i]
                    ril = rils[i]
                    data[ril][chrs][pos] = unit[i]
            else:
                unit = re.split(r'\t',line)
                unit[0] = 'Position'
                #print unit[0], unit[1], unit[2]
                rils.extend(unit)
   
    #for t in sorted(data['GN204']['Chr10'].keys(), key=int):
    #    print t
    return data

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
                mping = unit[0]
                data[mping] = [unit[1], unit[3]]
                #print mping, rils
    return data

# ""      "GN1"   "GN10"  "GN100" "GN101" "GN102" "GN103" "GN104" "GN105" "GN106" "GN107" "GN108
def read_rils(infile):
    rils = []
    data = defaultdict(lambda : defaultdict(lambda: defaultdict(lambda: str)))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            line = line.replace('"', '')
            if not line[0:1].isdigit():
                unit = re.split(r'\t',line)
                rils.extend(unit[1:])
    return rils

def subbam(mping, ril, bam, output):
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


def snp_idx_genotype(snp_idx, snpmap, ril, chrs):
    snp_gt = []
    for i in snp_idx:
        s = snpmap[ril][chrs][i]
        snp_gt.append(s)
    return snp_gt

def find_SNP_nearby(start, snpmap, ril, chrs):
    #array = [] 
    #array.extend(sorted(binmap[ril][chrs].keys(), key=int))
    print 'check', start, ril, chrs
    array_snp = []
    array_snp.extend(sorted(snpmap[ril][chrs].keys(), key=int))
    #mping after last bin on chromosome, return 0 mean genotype unknown
    try:
        if int(start) > int(array_snp[-1]):
            return ['NA']*10
    except:
        print 'check', start, ril, chrs
    #index = binarySearch(array, int(start))
    index_snp = binarySearch(array_snp, int(start))
    #print index_snp
    snp_idx   = array_snp[(index_snp[0]-10):(index_snp[0]+10)]
    snp_gt    = snp_idx_genotype(snp_idx, snpmap, ril, chrs)
    return snp_idx, snp_gt

    #check five flanking SNP if consistent
    #pos_gt = []
    #backward_snp_idx = array_snp[(index_snp[0]-4):(index_snp[0]+1)]
    #forward_snp_idx  = array_snp[index_snp[1]:(index_snp[0]+4)]
    #gt_1 = snp_type(backward_snp_idx, snpmap, ril, chrs)
    #gt_2 = snp_type(forward_snp_idx, snpmap, ril, chrs)


    #block genotype unknown, return
    #print 't:1 %s' %(index[1])
    #print 't: %s' %(binmap[ril][chrs][array[index[1]]])
    #if str(binmap[ril][chrs][array[index[1]]]) == 'NA':
    #    pos_gt.extend([array[index[1]], 'NA'])
    #else:
    #    #block genotype known, use snp to comfirm
    #    if gt_1[0] == gt_2[0] and str(gt_1[0]) == str(binmap[ril][chrs][array[index[1]]]):
    #        #snp consistent with genotype block
    #        pos_gt.extend([array[index[1]], gt_1[0]])
    #    else:
    #        #snp inconsistent with genotype block
    #        pos_gt.extend([array[index[1]], 'NA'])
    #return pos_gt

#NA      NA      0       0       0       0       0       0       0       0
def genotype_snp(gt):
    data = defaultdict(lambda : int())
    data_left5 = defaultdict(lambda : int())
    data_right5 = defaultdict(lambda : int())
    for i in range(0, len(gt)):
        data[gt[i]] += 1
        if i > 5 and i < 11:
            data_left5[gt[i]] += 1
        if i > 10 and i < 15:
            data_right5[gt[i]] += 1
    if data.has_key('1') and not data.has_key('0'):
    #pure HEG4 genotype
        if data_left5['1'] >= 2 and data_right5['1'] >= 2:
        # 5 SNP flanking both end also have HEG4 genotype
            return 1
        else:
            return '1'
    elif data.has_key('0') and not data.has_key('1'):
    #pure NB genotype
        if data_left5['0'] >= 2 and data_right5['0'] >= 2:
        # 5 SNP flanking both end also have NB genotype
            return 0
        else:
            return '0'
    elif data.has_key('0') and data.has_key('1'): 
    #recombination within block
        if data_left5.has_key('1') and not data_left5.has_key('0') and data_right5.has_key('1') and not data_right5.has_key('0'):
        # 5 SNP flanking both end have HEG4 genotype
            return 1
        elif data_left5.has_key('0') and not data_left5.has_key('1') and data_right5.has_key('0') and not data_right5.has_key('1'):
        # 5 SNP flanking both end have NB genotype
            return 0
        else:
            return '?'
    else:
        return '?'

def genotype_ping(mpings, rils, snpmap, outfile):
    ofile1 = open('%s.table.genotype.list' %(outfile), 'w')
    ofile2 = open('%s.table.ping_code.list' %(outfile), 'w')
    data = defaultdict(lambda : defaultdict(lambda : str()))
    for mping in sorted(mpings.keys()):
        print >> ofile1, '>%s' %(mping)
        #snp_header   = []
        #ril_genotype = []
        count = 0
        for ril in sorted(rils.keys(), key=int):
            count += 1
            #print '%s\t%s' %(mping, ril)
            ril_n = ril
            ril = 'GN%s' %(ril)
            
            p = re.compile(r'(\w+):(\d+)\-(\d+)')
            m = p.search(mping)
            chrs = ''
            start = 0
            end   = 0
            if m:
                chrs  = m.groups(0)[0]
                start = m.groups(0)[1]
                end   = m.groups(0)[2]
            #find_SNP_nearby(start, binmap, snpmap, ril, chrs)
            snp_idx, snp_gt = find_SNP_nearby(start, snpmap, ril, chrs)
            #snp_header = snp_idx
            #ril_genotype.append(snp_gt)
            ##output genotype and mPing
            if count == 1:
                snp_idx.insert(10, mping)
                print >> ofile1, '%s\t%s' %('SNP', '\t'.join(map(str, snp_idx)))
            snp_gt.insert(10, type)
            print >> ofile1, '%s\t%s' %(ril, '\t'.join(map(str, snp_gt)))
            ##genotype ping, turn NB genotype
            ping_gt = genotype_snp(snp_gt)
            if str(ping_gt) == '0':
                if mpings[mping][0] == 'NB':
                    ping_gt = '1'
            elif str(ping_gt) == '1':
                if mpings[mping][0] == 'NB':
                    ping_gt = '0'
 
            data[ril_n][mpings[mping][1]] = ping_gt
    ofile1.close()

    for ril in sorted(data.keys(), key=int):
        ping_codes1 = []
        ping_codes2 = []
        for c in sorted(data[ril].keys()):
            ping_codes1.append(data[ril][c])
            if str(data[ril][c]) == '1':
                ping_codes2.append(c)
        if len(ping_codes2) > 0:
            print >> ofile2, '%s\t%s\t%s\t%s\t%s' %(re.sub(r'GN', r'', ril), rils[re.sub(r'GN', r'', ril)], ''.join(map(str, ping_codes1)), str(len(ping_codes2)), ''.join(map(str, ping_codes2)))
        else:
            print >> ofile2, '%s\t%s\t%s\t0\tNA' %(re.sub(r'GN', r'', ril), rils[re.sub(r'GN', r'', ril)], ''.join(map(str, ping_codes1)))
    ofile2.close()

#return dict of ril->lib_name
def parse_bam_all(bam_list, r):
    data = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int())))
    for lib in sorted(re.split(r'\n', bam_list)):
        unit = re.split(r' |\t', lib)
        bam = os.path.split(unit[-1])[1]
        bam = re.sub(r'.recal.bam', r'', bam)
        bam = re.sub(r'.bam', r'', bam)
        #print lib, bam
        if r.search(bam):
            ril = r.search(bam).groups(0)[0]
            data[ril] = bam
            #print ril
    return data   

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--bin_map')
    parser.add_argument('--snp_map')
    parser.add_argument('--bam')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = 'ping.genotype'
    if not args.bin_map:
        args.bin_map = 'MPR.geno.bin'
    if not args.snp_map:
        args.snp_map = 'MPR.geno.data'
    if not args.bam:
        args.bam = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam'

    #bin map and snp genotype
    #binmap = convert_MAP(args.bin_map)
    snpmap = convert_MAP_SNP(args.snp_map)

    #read mping and rils
    pings = read_mping_list(args.input)
   
    #read ril
    #rils = read_rils(args.snp_map)
    bam_275 = subprocess.check_output('ls -all %s/*.bam' %(args.bam), shell=True)
    r2 = re.compile(r'RIL(\d+)\_') 
    ril_275 = parse_bam_all(bam_275, r2) 
    rils    = ril_275
    #genotyping 
    genotype_ping(pings, rils, snpmap, args.output)


if __name__ == '__main__':
    main()

