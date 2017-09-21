#!/opt/Python/2.7.3/bin/python
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
Python TraitDist.py --input ../input/test.trait

Calculate trait distance, which is 1-correlation for heatmap linakge cluster.

input:
Heading Days	Plant Height (cm) in Field	Biomass	Number of Tillers	Single Plant (Grain Yield) (g)	Height_Single_Plants	Grain_yield_per_tiller	Harvest_Index	seedling_12C	seedling_26C
103	92	120.2	22	42.0	95	1.88	0.299	0.8	0.87
99	99	206.7	30	90.0	100	3.03	0.259	0.857	0.915
93	94	198.9	32	86.7	97	2.68	0.413	NA	NA
96	109	219.6	29	87.8	112	3.06	0.384	0.71	0.915
93	90	202.1	33	65.2	94	1.96	0.438	0.777	0.9

output:
Heading Days	0.0	0.953779015647	1.00596447294	1.20958638286	1.06050052649	0.923654966779	0.911610554891	1.02818641826	0.829523583589	0.722077110754
Plant Height (cm) in Field	0.953779015647	0.0	0.553071785566	0.96100857909	0.684336320677	0.0431217117621	0.711427230417	0.924751011631	0.959409380786	0.799647853201
Biomass	1.00596447294	0.553071785566	0.0	0.331621474858	0.260184590036	0.529891430209	0.686228022061	0.929697127469	0.892338861216	0.737522257691
Number of Tillers	1.20958638286	0.96100857909	0.331621474858	0.0	0.662818841351	0.97396884445	1.29040822872	0.89747220721	1.03277812854	1.00796120936
Single Plant (Grain Yield) (g)	1.06050052649	0.684336320677	0.260184590036	0.662818841351	0.0	0.663066835402	0.212466308595	0.973428526742	0.831459468694	0.706865063153
Height_Single_Plants	0.923654966779	0.0431217117621	0.529891430209	0.97396884445	0.663066835402	0.0	0.679451688446	0.923504291385	0.934352740196	0.775346819125
Grain_yield_per_tiller	0.911610554891	0.711427230417	0.686228022061	1.29040822872	0.212466308595	0.679451688446	0.0	1.02822629063	0.811338434317	0.684037354541
Harvest_Index	1.02818641826	0.924751011631	0.929697127469	0.89747220721	0.973428526742	0.923504291385	1.02822629063	0.0	1.05567970349	1.01191933935
seedling_12C	0.829523583589	0.959409380786	0.892338861216	1.03277812854	0.831459468694	0.934352740196	0.811338434317	1.05567970349	0.0	0.560806546619
seedling_26C	0.722077110754	0.799647853201	0.737522257691	1.00796120936	0.706865063153	0.775346819125	0.684037354541	1.01191933935	0.560806546619	0.0


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
""      "GN1"   "GN10"  "GN100" "GN101" "GN102" "GN103" "GN104" "GN105" "GN106" "GN107" "GN108
"0100222046"    1       1       0       0       0       0       1       1       1       1
"0100500860"    1       1       0       0       0       0       1       1       1       1
'''
def dist(infile):
    data = defaultdict(list)
    dist1 = defaultdict(lambda : defaultdict(lambda : float))
    '''read bin table'''
    with open (infile, 'r') as filefh:
        header  = filefh.readline()
        header  = header.replace('"', '')
        headers = re.split('\t',header.rstrip())
        headers = headers[1:]
        #headers = sorted(headers)
        #print headers[0], headers[1]
        for line in filefh:
            line = line.rstrip()
            unit = re.split('\t',line)
            unit = unit[1:]
            for i in range(0,len(unit)):                
                value = re.sub(r'NA',r'NaN',unit[i])
                #print i, value 
                data[i].append(int(value))

    
    '''Pearson correlation''' 
    for rank in sorted(data.keys()):
        #print rank, ','.join(data[rank])
        for rank2 in sorted(data.keys()):
            #print rank, rank2
            list1, list2 = killnan(data[rank],data[rank2])
            dist0 = distance(list1,list2)
            #print '%s\t%s\t%s' % (rank, rank2, dist0)
            dist1[rank][rank2] = dist0
    '''output matirx, correlation and p-value'''
    #print "\t".join(headers)
    for rank in sorted(data.keys()):
        outline = []
        #outline.append(str(rank))
        outline.append(headers[rank])
        for rank2 in sorted(data.keys()):
            #print rank, rank2
            d = dist1[rank2][rank]
            outline.append(d)
        print "\t".join(map(str,outline)) 


def distance(x,y):
    p = len(x)
    m = sum(map(lambda (a,b): 1 if a == b else 0, zip(x,y)))
    return float(p-m)/p

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

