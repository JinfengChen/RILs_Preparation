#!~/BigData/software/miniconda2/bin/python
import sys
from collections import defaultdict
import re
import os
import argparse
from Bio import SeqIO
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.spatial.distance as distance
import matplotlib as mp
mp.use('pdf')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pandas as pd


def usage():
    test="name"
    message='''
python Trait_Tree_Label.py --input ../input/May28_2013.RIL.trait.table.QTL.trait.txt.ALL.dist

    '''
    print message

def set_ticks_XY(ax, ypos, ylim, chrs):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(True)
    ax.yaxis.set_ticks_position('right')
    ax.xaxis.set_ticks(chrs[2], minor=True)

    # put the major ticks at the middle of each cell
    ax.set_yticks(ypos, minor=False)
    ax.set_xticks(chrs[1], minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels(chrs[0], minor=False)
    ax.set_yticklabels(ylim, minor=False)

    # rotate the
    #plt.xticks(rotation=0)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()


    for t in ax.xaxis.get_major_ticks():
    #    t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
 

    return ax


def set_ticks_XY_Right(ax, ypos, ylim):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(True)
    ax.yaxis.set_ticks_position('right')

    # put the major ticks at the middle of each cell
    ax.set_yticks(ypos, minor=False)
    #ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels([], minor=False)
    ax.set_yticklabels(ylim, minor=False)

    # rotate the
    #plt.xticks(rotation=0)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax


def set_ticks_XY_empty0(ax):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    #ax.set_yticks(np.arange(xlen) + 0.5, minor=False)
    #ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels([], minor=False)
    ax.set_yticklabels([], minor=False)

    # rotate the
    plt.xticks(rotation=90)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax


def set_ticks_XY_empty(ax):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    #ax.set_yticks(np.arange(xlen) + 0.5, minor=False)
    #ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels([], minor=False)
    #ax.set_yticklabels([], minor=False)

    # Set y ticklabels font size
    for tck in ax.get_yticklabels():
        tck.set_fontsize(2)

    # rotate the
    plt.xticks(rotation=90)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax



def LOD_curve(mat, labs, bins):
    # clustering
    dist_mat = mat
    linkage_matrix = linkage(dist_mat, "single")

    # Create a figure.
    figsize=(20,10)
    fig = plt.figure(figsize=figsize)
    plt.subplots_adjust(bottom=0.2)
  
    # split plot into 4 rows X 2 cols. gs[1, 0] is 2th row and 1th column, which is left bottom region.
    # gs[0:, 0]  is all rows of 1th column
    ncurve = 272.00
    gs = gridspec.GridSpec(int(ncurve), 2, wspace=0.04, hspace=0.15, width_ratios=[0.15,1], height_ratios=[1/ncurve]*int(ncurve))

    #dendrogram
    ax0 = fig.add_subplot(gs[0:, 0])
    ddata = dendrogram(linkage_matrix,
                   color_threshold=10,
                   orientation='left',
                   labels=labs)

    ax0 = set_ticks_XY_empty(ax0)
    fig.savefig('%s.pdf' %('Trait_Tree_Label'), bbox_inches='tight')
'''
    ymax = 10
    ymin = 0
    #xmax = 373245519
    xmax = int(max(bins['Position']))
    xmin = 0 
    count = 0
    for i in ddata['leaves']:
        #print '%s\t%s' %(i, labs[i])
        ax1 = fig.add_subplot(gs[int(ncurve)-count-1, 1])
        positions = bins['Position']
        genotypes = bins[labs[i]]
        lastpos = 0
        lastgt  = 2
        start   = 0
        x1range = []
        y1range = (0, 10)
        c1range = []
        ax1.set_xlim([xmin, xmax])
        ax1.set_ylim([ymin, ymax])
        for j in range(len(positions)):
            if (lastgt == 2):
                lastgt = int(genotypes[j])
                lastpos = positions[j]
            elif(int(genotypes[j]) != lastgt):
                r = (start, lastpos)
                start = lastpos+1
                lastgt = int(genotypes[j])
                lastpos = positions[j]
                x1range.append(r)
                color = 'blue' if int(lastgt) == 1 else 'red'
                c1range.append(color)
            else:
                lastpos = positions[j]
        print i, labs[i], x1range
        ax1.broken_barh(x1range, y1range, facecolors=tuple(c1range), edgecolor = "none")
        ax1=set_ticks_XY_empty0(ax1) 
        count += 1
    # set axis
    #ax0 = set_ticks_XY_empty(ax0)
    fig.savefig('%s.pdf' %('QTL_BinMap_Label'), bbox_inches='tight')'''


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


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

    #lodfile = '../input/MPR.cross.uniq.QTL.mr.table.new'
    #lodcutfile = '../input/MPR.cross.uniq.QTL.mr.table.LOD_threshold.new'
    #chrfile = '../input/MSU7.Chr.midpoint'
    pdist = pd.read_table(args.input, index_col=0, header=None)
    #midpoint = pd.read_table(chrfile, header=None)
    binmap   = pd.read_table('MPR.geno.bin.uniq.new')
    #lodcut= pd.read_table(lodcutfile)   
 
    mat = np.array([[100,  80,  50],
                [1,  0.8, 0.5],
                [30, 80,  80],
                [0.5,  0.8, 1]])

    #pairwise_dists = distance.squareform(distance.pdist(mat))
    LOD_curve(pdist, pdist.index, binmap)

if __name__ == '__main__':
    main()

