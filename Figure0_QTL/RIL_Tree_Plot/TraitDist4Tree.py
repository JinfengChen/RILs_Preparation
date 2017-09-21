#!~/BigData/software/miniconda2/bin/python
import sys
from collections import defaultdict
from numpy import *
from scipy.stats.stats import pearsonr
import re
import os
import argparse

def usage():
    test="name"
    message='''
python TraitDist4Tree.py --input ../input/May28_2013.RIL.trait.table.QTL.trait.txt.ALL > ../input/May28_2013.RIL.trait.table.QTL.trait.txt.ALL.dist

input:
Sample  Heading Days    Plant Height (cm) in Field      Biomass Number of Tillers       Single Plant (Grain Yield) (g)  Non-ref mPing   Unique mPing    Ping    Cold_ratio_index        CV(cold_ratio_index)
GN-1    103     92      120.2   22      42.0    200     25      5       0.92    12.50   243     94      86      95      1.88    0.299   0.8     0.87    
GN-2    99      99      206.7   30      90.0    128     11      2       0.935   13.68   269     91      81      100     3.03    0.259   0.857   0.915   
GN-3    93      94      198.9   32      86.7    252     13      6       NA      NA      NA      87      89      97      2.68    0.413   NA      NA      
GN-4    96      109     219.6   29      87.8    198     26      5       0.777   29.77   149     90      99      112     3.06    0.384   0.71    0.915

output:
GN1     0.0     0.480446927374  0.550465549348  0.489757914339  0.481936685289  0.462569832402  0.467039106145  0.458100558659  0.460335195531  0.395903165736  0.494227188082  0.633519553073  0.435754189944  0.55232774
GN10    0.480446927374  0.0     0.452886405959  0.554562383613  0.505027932961  0.522905027933  0.431284916201  0.525139664804  0.488640595903  0.444320297952  0.509124767225  0.533705772812  0.443202979516  0.500186219739
GN100   0.550465549348  0.452886405959  0.0     0.450279329609  0.553445065177  0.506517690875  0.566108007449  0.496834264432  0.500558659218  0.461452513966  0.439851024209  0.503165735568  0.511731843575  0.56350093109
GN101   0.489757914339  0.554562383613  0.450279329609  0.0     0.522532588454  0.535195530726  0.553817504655  0.465176908752  0.459962756052  0.548975791434  0.537802607076  0.534823091248  0.512849162011  0.4908752327
GN102   0.481936685289  0.505027932961  0.553445065177  0.522532588454  0.0     0.574301675978  0.342644320298  0.483426443203  0.433519553073  0.521787709497  0.55009310987   0.389944134078  0.56834264432   0.642458100

    '''
    print message

def killnan(list1, list2):
    list3 = []
    list4 = []
    for i in range(len(list1)):
        if isnan(list1[i]) or isnan(list2[i]):
            continue
        else:
            #print list1[i], list2[i]
            list3.append(list1[i]) 
            list4.append(list2[i]) 
    #print len(list1), len(list2), len(list3), len(list4)
    return list3, list4 


'''
Sample  Heading Days    Plant Height (cm) in Field      Biomass Number of Tillers       Single Plant (Grain Yield) (g)  Non-ref mPing   Unique mPing    Ping    Cold_ratio_index        CV(cold_ratio_index)
GN-1    103     92      120.2   22      42.0    200     25      5       0.92    12.50   243     94      86      95      1.88    0.299   0.8     0.87    
GN-2    99      99      206.7   30      90.0    128     11      2       0.935   13.68   269     91      81      100     3.03    0.259   0.857   0.915   
'''
def dist(infile):
    data = defaultdict(list)
    dist1 = defaultdict(lambda : defaultdict(lambda : float))
    '''read trait table'''
    with open (infile, 'r') as filefh:
        header = filefh.readline()
        #print headers[0], headers[1]
        for line in filefh:
            line = line.rstrip()
            unit = re.split('\t',line)
            ril  = unit[0]
            unit = unit[1:]
            for i in range(0,len(unit)):                
                value = re.sub(r'NA',r'NaN',unit[i])
                #print ril, i, value 
                data[ril].append(float(value))

    
    '''Pearson correlation''' 
    for ril1 in sorted(data.keys()):
        #print rank, ','.join(data[rank])
        for ril2 in sorted(data.keys()):
            #print rank, rank2
            list1, list2 = killnan(data[ril1],data[ril2])
            dist0 = distance(list1,list2)
            #print '%s\t%s\t%s' % (ril1, ril2, dist0)
            dist1[ril1][ril2] = dist0
    '''output matirx, correlation and p-value'''
    #print "\t".join(headers)
    for ril1 in sorted(data.keys()):
        outline = []
        #outline.append(str(rank))
        outline.append(ril1)
        for ril2 in sorted(data.keys()):
            #print rank, rank2
            d = dist1[ril2][ril1]
            outline.append(d)
        print "\t".join(map(str,outline)) 


def distance(x,y):
    p = pearsonr(x, y)
    d = 1-float(p[0])
    return d
    #p = len(x)
    #m = sum(map(lambda (a,b): 1 if a == b else 0, zip(x,y)))
    #return float(p-m)/p

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    dist(args.input)

if __name__ == '__main__':
    main()

