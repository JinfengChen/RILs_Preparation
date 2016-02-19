#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random

def usage():
    test="name"
    message='''
python ../../bin/Pseudo_TEinsertion_Genome.py --repeat mPing_Ping_Pong.fa --gff HEG4.ALL.mping.non-ref.gff --genome MSU_r7.fa

--repeat: fasta file contains sequence of TEs in gff file (no need to have TSD=... after sequence name).
--gff:    gff file contains location of non_reference TE insertions in reference genome. TE=name;TSD=TAA; must present in gff, so we know which TE to insert and the length of TSD.
--genome: reference genome sequence
--project: prefix of output files

The script make a pseudogenome by inserting non_reference TEs into reference chromosome. A new gff will be generated which contains the locations of TEs in pseudogenome. The new pseudogenome can be used to analyzing DNA methylation and histone modification of TEs and flanking regions.
    '''
    print message

##
def write_fasta(ref, filename):
    ofile = open(filename, "w")
    for chrn in sorted(ref.keys()):
        seq = Seq(ref[chrn])
        newrecord = SeqRecord(seq, id=chrn,description="")
        SeqIO.write(newrecord, ofile, "fasta")
    ofile.close()

##repeat fasta
def fasta_te(fastafile, te):
    fastaid = defaultdict(str)
    fastatsd = defaultdict(list)
    ids = []
    for record in SeqIO.parse(fastafile, "fasta"):
        fastaid[record.id] = str(record.seq)
        ids.append(str(record.id))
    #get tsd
    s = re.compile(r'>(.\S+)\s+TSD\=(.*)')
    disc = SeqIO.index(fastafile, "fasta")
    for d in ids:
        lines = re.split('\n', disc.get_raw(d).decode())
        head  = lines[0]
        m = s.search(head)
        repid = ''
        tsd   = ''
        if m:
            repid = m.groups(0)[0]
            tsd   = m.groups(0)[1]
            fastatsd[repid] = [fastaid[repid], tsd]
            #print repid, tsd
    return fastatsd

##reference fasta
def fasta_ref(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile, "fasta"):
        fastaid[record.id] = str(record.seq)
        print str(record.id), str(len(str(record.seq)))
    return fastaid

##random choose pos from chromosome
def random_pos(chrn, ref, repid, repseq, reptsd):
    seq_len = len(ref[chrn])
    rand_pos= random.randint(1,seq_len)
    strand  = 1 if random.randint(1,10) > 5 else 0
    #tsd
    tsdstart = rand_pos - len(reptsd)
    chrseq   = ref[chrn]
    tsdseq   = chrseq[tsdstart:rand_pos]
    #tsdseq   = chrseq[tsdstart:rand_pos] if strand == 1 else str(Seq(chrseq[tsdstart:rand_pos]).reverse_complement())
    return [chrn, rand_pos, rand_pos+1, strand, repid, repseq, tsdseq]

##except for the first insertion, all the following insertion need to add inserted sequence to their position
def update_pos(data):
    ##tsd only add one copy of itself
    #add_len = len(repseq) + 1*len(tsd)
    add_total = defaultdict(int)
    for i in range(0, len(data)):
        chrn  = data[i][0]
        start = data[i][1]
        end   = data[i][2]
        repid = data[i][4] 
        repseq= data[i][5]
        reptsd= data[i][6]
        add_len = len(repseq) + 1*len(reptsd)
        
        #print chrn, start, end, repid, repseq, reptsd
        if add_total.has_key(chrn):
            # add previous insertion length to the coordinate
            data[i][1] = start + add_total[chrn]
            data[i][2] = end   + add_total[chrn]
            add_total[chrn] += add_len
        else:
            # first insertion, no change to the coordinate
            add_total[chrn] = add_len
        #print data[i][0], data[i][1], data[i][2]
    return data

##insertion element into genome
def insert_element(pos, ref, ofile):
    chrn = pos[0]
    chrseq = ref[chrn]
    half1 = chrseq[:pos[2]]
    half2 = chrseq[pos[2]:]
    repid =pos[4]
    repseq=pos[5]
    reptsd=pos[6]
    repname=pos[7]
    repend =int(pos[1]) + len(repseq) + 1*len(reptsd) - 1
    #print 'insert: %s, %s' %(reptsd, repseq)
    #Chr1    not.give        transposable_element_attribute  1132975 1132977 -       .       .       ID=Chr1.1132977.spanners;avg_flankers=17;spanners=0;type=homozygous;TE=mping;TSD=TAA
    gff_newline = '%s\tPseudoGenome\tTransposable_element\t%s\t%s\t%s\t.\t.\tID=%s_%s_%s;Original_ID=%s;TE=%s;TSD=%s;' %(chrn, pos[1]+len(reptsd), repend, pos[3], chrn, pos[1]+len(reptsd), repend, repid, repname, reptsd)
    print >> ofile, gff_newline
    ##we choose sequence at target site as tsd, not use tsd provided
    tsdstart = pos[2] - len(reptsd)
    tsdseq   = chrseq[tsdstart:pos[2]]
    newseq   = ''
    if pos[3] == '+':
        newseq = half1 + repseq + tsdseq + half2
        #print tsdseq, repseq
    else:
        repseq_seq = Seq(repseq)
        repseq_rec = repseq_seq.reverse_complement()
        #print tsdseq, str(repseq_rec)
        newseq = half1 + str(repseq_rec) + tsdseq + half2
    ref[chrn] = newseq

def gff_parser(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd: 
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('#'): 
                unit  = re.split(r'\t',line)
                start = int(unit[3]) 
                end   = int(unit[4])
                chro  = unit[0]
                strand= unit[5]
                temp  = defaultdict(str)
                attrs = re.split(r';', unit[8])
                for attr in attrs:
                    #print attr
                    if not attr == '':
                        #print 'yes'
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value
                repid   = temp['ID'] if temp.has_key('ID') else '%s_%s_%s' %(chro, start, end)
                repname = temp['TE']
                reptsd  = temp['TSD'] if temp.has_key('TSD') else 'TAA'
                #repfam  = temp['Class']
                data[chro].append([start, end, strand, repid, repname, reptsd, unit[8]])
                #print '%s\t%s\t%s\t%s\t%s\t%s' %(repid, chro, start, end, repname, reptsd)
                #print '%s\t%s\t%s' %(chro, start, end)
    return data

##write simulated insertion into gff format
def writegff(data, gff_out):
    ofile = open(gff_out, 'w') 
    for pos in data:
        chrn   = pos[0]
        start  = pos[1]
        end    = pos[2]
        strand = '+' if pos[3] == 1 else '-'
        repid  = pos[4]
        repseq = pos[5]
        reptsd = pos[6]
        te_id  = '%s.%s.%s' %(repid, chrn, str(start))
        print >> ofile, '%s\tMSU7\t%s\t%s\t%s\t.\t%s\t.\tID=%s;TSD=%s;' %(chrn, repid, str(start), str(end), strand, te_id, reptsd)

def pick_te(element, te):
    if te == 'ALL':
        te_ids = sorted(element.keys())
        te_num = len(element.keys())
        #print te_num
        index  = random.randint(1, int(te_num))
        te_inf = [te_ids[index-1], element[te_ids[index-1]][0], element[te_ids[index-1]][1]]
    else:
        te_inf = [te, element[te][0], element[te][1]]
    return te_inf

##main function of simulation
def simulate(ref, element, gff, prefix):
    data = []
    for chro in sorted(gff.keys()): 
        for ins in sorted(gff[chro]):
            #print chro, ins[0], ins[1], ins[2], ins[3], ins[4], element[ins[4]], ins[5]
            data.append([chro, ins[0], ins[1], ins[2], ins[3], element[ins[4]], ins[5], ins[4]])
    data = sorted(data, key = lambda x: (x[0], x[1]))

    ###update insertion to new pseudogenome position
    gff_out    = '%s.gff' %(prefix)
    genome_out = '%s.fa' %(prefix)
    datac      = data
    datac      = update_pos(datac) 
    
    ###insert TE into reference genome
    ofile = open(gff_out, 'w')
    for pos in datac:
        insert_element(pos, ref, ofile)
    ofile.close()
    write_fasta(ref, genome_out)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--repeat')
    parser.add_argument('--genome')
    parser.add_argument('--gff')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.repeat) > 0 or len(args.genome) > 0 or len(args.gff) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = '%s.Pseudo_mPing' %(os.path.split(os.path.splitext(args.genome)[0])[1])
    pseudogenome = '%s.fa' %(args.project)
    pseudogff    = '%s.gff' %(args.project) 
    print pseudogenome
    print pseudogff

    ref    = fasta_ref(args.genome)
    repeat = fasta_ref(args.repeat)
    gff    = gff_parser(args.gff)
    simulate(ref, repeat, gff, args.project)

if __name__ == '__main__':
    main()

