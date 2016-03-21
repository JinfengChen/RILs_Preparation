import sys
import HTSeq
import numpy
import matplotlib as mpl
mpl.use('pdf')
from matplotlib import pyplot

#sortedbamfile = HTSeq.BAM_Reader( "../input/DHS.Chr1.unique.bam" )
#sortedbamfile = HTSeq.BAM_Reader( "../input/DHS.unique.bam" )
#gtffile = HTSeq.GFF_Reader( "../input/MSU7.gene.exon_number.gtf" )
sortedbamfile = HTSeq.BAM_Reader(sys.argv[1])
gtffile = HTSeq.GFF_Reader(sys.argv[2])

halfwinwidth = 2000
fragmentsize = 150
#total = 60745783.00/1000000 ## nucleosome
#total = 7480914/1000000 ## nucleosome chr1
#total = 23299296/1000000 #DHS unique
#gsize = 372000000

#coverage = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )
#for almnt in bamfile:
#   if almnt.aligned:
#      #almnt.iv.length = fragmentsize
#      print almnt.iv
#      if not almnt.iv.start < 500:
#          coverage[ almnt.iv ] += 1

tsspos = set()
for feature in gtffile:
   if feature.type == "exon" and feature.attr["exon_number"] == "1":
      #print feature.iv.start_d_as_pos.pos
      if feature.iv.start_d_as_pos.pos > 5000:
          tsspos.add( feature.iv.start_d_as_pos )

profile = numpy.zeros( 2*halfwinwidth, dtype='i' )      
for p in tsspos:
   window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth - fragmentsize, p.pos + halfwinwidth + fragmentsize, "." )
   for almnt in sortedbamfile[window]:
       almnt.iv.length = fragmentsize
       if p.strand == "+":
           start_in_window = almnt.iv.start - p.pos + halfwinwidth 
           end_in_window   = almnt.iv.end   - p.pos + halfwinwidth 
       else:
           start_in_window = p.pos + halfwinwidth - almnt.iv.end
           end_in_window   = p.pos + halfwinwidth - almnt.iv.start
       start_in_window = max( start_in_window, 0 )
       end_in_window = min( end_in_window, 2*halfwinwidth )
       if start_in_window >= 2*halfwinwidth or end_in_window < 0:
           continue
       profile[ start_in_window : end_in_window ] += 1

ofile = open('DHS_TSS.profile', 'w')
for i in range(len(profile)):
    print >> ofile, profile[i]
ofile.close()

pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )  
pyplot.savefig('DHS_TSS.pdf')
#pyplot.savefig('%s.pdf' %(sys.argv[2]))

