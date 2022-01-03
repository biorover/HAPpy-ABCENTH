#!/usr/bin/env python

import argparse,sys,os
from intervaltree import IntervalTree
import pandas as pd

def chose_clusters(hittabm,maxintron = 10000):
    """
    reads through a hit table to find best-scoring non-overlapping clusters of hits
    """
    hittab.columns = ['cluster','seqid','pid','dunno','dontcare','strand','qstart','qend','sstart','send',
                  'evalue','score','pstart','pend','qlen','ei','en']
    hittab = hittab.sort_values(['seqid','sstart']) 
    clusters = list(set(hittab['cluster']))
    seqs = list(set(hittab['seqid']))
    seqdict = {seq:IntervalTree() for seq in seqs}
    loci = []
    cols = ['cluster','seqid','start','end','score']
    locrowis = {}
    for cluster in clusters:
        cdf = hittab.query('cluster == "' + cluster + '"')
        strand,ei,seqid,lastend = '.',0,'.',0
        for i,row in cdf.iterrows():
            start = min((row['sstart'],row['send']))
            end = max((row['sstart'],row['send']))        
            if row['seqid'] == seqid and row['strand'] == strand and end - lastend < maxintron and \
                    ((strand == "+" and row['ei'] > ei) or (strand == '-' and row['ei'] < ei)):
                loci[-1][3] = end
                loci[-1][4] += row['score']
                locrowis[len(loci)].append(i)
            else:
                loci.append([row['cluster'],row['seqid'],start,end,row['score']])
                locrowis[len(loci)] = [i]
            strand,ei,seqid,lastend = row['strand'],row['ei'],row['seqid'],end

    locdf = pd.DataFrame(loci)
    locdf.columns = cols
    locdf = locdf.sort_values('score',ascending = False)
    keepis = []
    for index,row in locdf.iterrows():
        start = min((row['start'],row['end']))
        end = max((row['start'],row['end']))
        if not seqdict[row['seqid']].overlaps(start,end):
            seqdict[row['seqid']][start:end] = 1
            keepis.extend(locrowis[index])
    outdf = hittab.loc[keepis,:]
    outdf = outdf.sort_values(['seqid','sstart'])
    return outdf

def main():
    parser = argparse.ArgumentParser("ClusterChooser: reads through a hit table to find best-scoring non-overlapping clusters of hits")
    parser.add_argument('--table', help = 'abcenth format hit table')
    parser.add_argument('--max_intron', help = 'maximum intron size (default = 10000', default = 10000, type = int)
    args = parser.parse_args()

    hittab = pd.read_csv(args.table,sep="\t",header=None)
    hittab = chose_clusters(hittab,args.max_intron)
    hittab.to_csv('/dev/stdout', sep = '\t', header = False, index = False)

if __name__ == "__main__":
    main()