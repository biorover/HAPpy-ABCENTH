#!/usr/bin/env python

import argparse,sys,os
from intervaltree import IntervalTree
import pandas as pd

def generate_clusterdict(cluster_info_df):
    """
    reads through cluster info and generates a dictionary of exon phases for each cluster
    """
    cdict = {}
    cluster_info_df.columns = ['cluster','exon','exonN','sphase','ephase','qlen']
    cdf = cluster_info_df.sort_values(['cluster','exon'])
    for index,row in cdf.iterrows():
        if not row['cluster'] in cdict:
            cdict[row['cluster']] = ""
        cdict[row['cluster']] += str(row[sphase])
    return cdict

def chose_clusters(hittab,cluster_dict,maxintron = 10000, criterion = 'matches'):
    """
    reads through a hit table to find best-scoring non-overlapping clusters of hits
    """
    hittab.columns = ['cluster','seqid','pid','dunno','dontcare','strand','qstart','qend','sstart','send',
                  'evalue','score','pstart','pend','qlen','ei','en']
    
    hittab = hittab.sort_values(['seqid','sstart'])
    hittab['exonstructure'] = hittab['cluster'].map(cluster_dict)
    estructs = list(set(hittab['exonstructure']))
    seqs = list(set(hittab['seqid']))
    seqdict = {seq:IntervalTree() for seq in seqs}
    loci = []
    cols = ['cluster','seqid','start','end','score','matches']
    locrowis = {}
    for estruct in estructs:
        cdf = hittab.query('exonstructure == "' + estruct + '"')
        strand,ei,seqid,lastend = '.',0,'.',0
        for i,row in cdf.iterrows():
            start = min((row['sstart'],row['send']))
            end = max((row['sstart'],row['send']))
            idb = row['qlen'] * row['pid']

            if row['seqid'] == seqid and row['strand'] == strand and end - lastend < maxintron and \
                    ((strand == "+" and row['ei'] > ei) or (strand == '-' and row['ei'] < ei)):
                loci[-1][3] = end
                loci[-1][4] += row['score']
                loci[-1][5] += idb
                locrowis[len(loci)].append(i)
            elif row['seqid'] == seqid and row['strand'] == strand and end - lastend < maxintron and \
                    row['ei'] == ei:
                if lastscore < row['score']:
                    loci[-1][4] += row['score'] - lastscore
                    loci[-1][3] = end
                    locrowis[len(loci)][-1] == i
            else:
                loci.append([row['cluster'],row['seqid'],start,end,row['score'],idb])
                locrowis[len(loci)] = [i]
            strand,ei,seqid,lastend,lastscore = row['strand'],row['ei'],row['seqid'],end,row['score']

    locdf = pd.DataFrame(loci)
    locdf.columns = cols
    locdf = locdf.sort_values(criterion,ascending = False)
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
    parser.add_argument('--cluster_info', help = 'cluster info from ABCENTH.HitTabFilter')
    parser.add_argument('--criterion', help = 'criterion for cluster selection (can be matches (default) or score)', default = 'matches')
    args = parser.parse_args()

    hittab = pd.read_csv(args.table,sep="\t",header=None)
    clustertab = pd.read_csv(args.cluster_info,sep='\t',header=None)
    clusterdict = generate_cluster_dict(clustertab)
    hittab = chose_clusters(hittab,clusterdict,args.max_intron,args.criterion)
    hittab.to_csv('/dev/stdout', sep = '\t', header = False, index = False)

if __name__ == "__main__":
    main()