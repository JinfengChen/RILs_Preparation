import sys
import HTSeq
import numpy
import matplotlib as mpl
mpl.use('pdf')
from matplotlib import pyplot


def get_nucleosome_profile(tsspos, upwinwidth, downwinwidth, coverage):
    profile = numpy.zeros( 2*halfwinwidth, dtype=float )
    for p in tsspos:
        window = HTSeq.GenomicInterval( p.chrom, p.pos - upwinwidth, p.pos + downwinwidth, "." )
        wincvg = numpy.fromiter( coverage[window], dtype='i', count=upwinwidth + downwinwidth)
        if p.strand == "+":
            profile += wincvg
        else:
            profile += wincvg[::-1]
    return profile

def write_profile(profile, outfile):
    ofile = open(outfile, 'w')
    for i in range(len(profile)):
        profile[i] = float(profile[i])
        print >> ofile, profile[i]
    ofile.close()
    

def plot_profile(profile, upwinwidth, downwinwidth, outfile):
    pyplot.figure()
    pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile)
    #pyplot.ylim((0,0.00002)) 
    pyplot.savefig(outfile)


bamfile = HTSeq.BAM_Reader(sys.argv[1])
gtffile = HTSeq.GFF_Reader(sys.argv[2])
#bamfile = HTSeq.BAM_Reader( "../input/Nucleosome.Chr1.unique.bam" )
#gtffile = HTSeq.GFF_Reader( "../input/MSU7.gene.exon_number.HighExp.Chr1.gtf" )

halfwinwidth = 2000
fragmentsize = 73
total = 60745783.00/1000000 ## nucleosome
gsize = 20478.00*2000.00
readlen = 36

coverage = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )
for almnt in bamfile:
   if almnt.aligned:
      #almnt.iv.length = fragmentsize
      if almnt.iv.strand == '+':
          almnt.iv.start = almnt.iv.start + fragmentsize - 10
          almnt.iv.end = almnt.iv.start + 20
#          if not almnt.iv.start < 500:
#              coverage[ almnt.iv ] += 1
      else:
          almnt.iv.start = almnt.iv.start + readlen - fragmentsize - 10
          almnt.iv.end = almnt.iv.start + 20
      print almnt.iv
      if not almnt.iv.start < 500:
          coverage[ almnt.iv ] += 1

tsspos = set()
ttspos = set()
for feature in gtffile:
   if feature.type == "exon" and feature.attr["exon_number"] == "1":
      #print feature.iv.start_d_as_pos.pos
      if feature.iv.start_d_as_pos.pos > 5000:
          tsspos.add( feature.iv.start_d_as_pos )
   elif feature.type == "exon" and feature.attr["exon_number"] == "-1":
      if feature.iv.start_d_as_pos.pos > 5000: 
          ttspos.add( feature.iv.end_d_as_pos )

'''
exon_s = set()
exon_e = set()
for feature in gtffile:
   if feature.type == "exon":
       if feature.iv.start_d_as_pos.pos > 5000:
           if feature.iv.strand == '+':
               #print feature.iv.start_d_as_pos
               exon_s.add(feature.iv.start_d_as_pos)
               exon_e.add(feature.iv.end_d_as_pos)
           else:
               exon_s.add(feature.iv.end_d_as_pos)
               exon_e.add(feature.iv.start_d_as_pos)
'''

#profile_tss = numpy.zeros( 2*halfwinwidth, dtype=float )      
#for p in tsspos:
#   window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
#   wincvg = numpy.fromiter( coverage[window], dtype='i', count=2*halfwinwidth )
#   if p.strand == "+":
#      profile_tss += wincvg
#   else:
#      profile_tss += wincvg[::-1]

##TSS profile
profile_tss = get_nucleosome_profile(tsspos, halfwinwidth, halfwinwidth, coverage)
write_profile(profile_tss, 'Nucleosome_TSS.profile')
plot_profile(profile_tss, halfwinwidth, halfwinwidth, 'Nucleosome_TSS.pdf')

##TTS_profile
profile_tts = get_nucleosome_profile(ttspos, halfwinwidth, halfwinwidth, coverage)
write_profile(profile_tts, 'Nucleosome_TTS.profile')
plot_profile(profile_tts, halfwinwidth, halfwinwidth, 'Nucleosome_TTS.pdf')

##Exon_profile
#profile_exon_s = get_nucleosome_profile(exon_s, 2000, 2000, coverage)
#write_profile(profile_exon_s, 'Nucleosome_Exon_S.profile')
#plot_profile(profile_exon_s, 2000, 2000, 'Nucleosome_Exon_S.pdf')

#profile_exon_e = get_nucleosome_profile(exon_e, 2000, 2000, coverage)
#write_profile(profile_exon_e, 'Nucleosome_Exon_E.profile')
#plot_profile(profile_exon_e, 2000, 2000, 'Nucleosome_Exon_E.pdf')


#ofile = open('Nucleosome.profile', 'w')
#for i in range(len(profile)):
#    profile[i] = float(profile[i])/(gsize*total)
#    print >> ofile, profile[i]
#ofile.close()

#pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile/total )
##pyplot.ylim((0,0.00002)) 
#pyplot.savefig('example01_chr1_20bp_bothend.pdf')

