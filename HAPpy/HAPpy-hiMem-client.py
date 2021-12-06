#!/usr/bin/env python

import sys
import os
import shutil
import argparse
from HAPpy import genome_fork as genome
from HAPpy import toolbox_for_HAP
import subprocess
import ete3
import time
import copy
import threading
import shlex
import tempfile
from collections import defaultdict,deque
from intervaltree import IntervalTree,Interval
import pickle

def output_fastas(clusters, prot_seq_dict, search_mode, fasta_dir,ref_genome_list = None,annotations = None):
    """Uses clusters dictionary built by "build_clusters" to break input sequences into fastas\
    for individual clusters. If search_mode = "exons", builds fastas for each exon."""
    sys.stderr.write('writing fastas for ' + str(len(clusters)) + ' clusters\n')
    cluster_dict = {}
    lengths_dict = {}
    if search_mode == 'exons':
        ref_genomes = []
        for i in range(len(ref_genome_list)):
            ref_genomes.append(genome.Genome(ref_genome_list[i]))
            ref_genomes[-1].read_gff(annotations[i])
    for cluster_index in range(len(clusters)):
        cluster = clusters[cluster_index]
        exon_file_root = fasta_dir + "/cluster" + str(cluster_index)
        cluster_dict[cluster_index] = []
        exon_number = None
        exon_phases = None
        exon_lengths = None
        for seqname in cluster:
            transcript_id = "_".join(seqname.split("_")[1:]) # because mafft adds "{number}_" to the front of the seq names
            if search_mode == "fl":
                seq_out = open(exon_file_root + '.fa','a')
                seq_out.write(">" + transcript_id + "\n" + prot_seq_dict[transcript_id] + '\n')
                seq_out.close()
            elif search_mode == "exons":
                for ref_genome in ref_genomes:
                    try:
                        cluster_dict[cluster_index].append(ref_genome.annotations.transcript[transcript_id])
                    except KeyError:
                        continue
                    ##This doesn't really fit the new program logic but might get re-incorporated depending on how I
                    ##decide to handle exon-based candidate loci parsing
                    #seq_out = open(exon_file_root + 'fullLenForHMM.fa','a')
                    #seq_out.write(ref_genome.annotations.transcript[transcript_id].get_fasta(seq_type='protein') + '\n')
                    #seq_out.close()
                    if not exon_number:
                        exon_number = len(ref_genome.annotations.transcript[transcript_id].child_list)
                    else:
                        if len(ref_genome.annotations.transcript[transcript_id].child_list) != exon_number:
                            sys.stderr.write("uh oh, transcripts with varying number of exons in: cluster " + str(cluster_index) + ", specifically gene " + transcript_id)
                            continue
                    sortlist = []
                    for CDS_id in ref_genome.annotations.transcript[transcript_id].child_list:
                        sortlist.append((ref_genome.annotations.CDS[CDS_id].get_coords(),CDS_id))
                    sortlist.sort()
                    if ref_genome.annotations.transcript[transcript_id].strand == "-":
                        sortlist.reverse()
                    if not exon_phases:
                        exon_phases = []
                        exon_lengths = []
                        current_phase = [0,0]
                        for cds_info in sortlist:
                            current_phase[1] = (current_phase[0] + cds_info[0][1] - cds_info[0][0] + 1) % 3
                            exon_lengths.append((cds_info[0][1] - current_phase[0] - cds_info[0][0] + 1) / 3)
                            exon_phases.append(current_phase[:])
                            current_phase[0] = current_phase[1]
                    for CDS_index in range(len(sortlist)):
                        cds_out = open(exon_file_root + "exon" + str(CDS_index) + 'of' + str(exon_number - 1) + 'phases' + str(exon_phases[CDS_index][0]) + 'and' + str(exon_phases[CDS_index][1])+ '.fa','a')
                        CDS_id = sortlist[CDS_index][1]
                        phase_slice_offset = (3 - exon_phases[CDS_index][0]) % 3
                        cds_out.write(">" + transcript_id + "-" + str(CDS_index) + '\n' + genome.Sequence(ref_genome.annotations.CDS[CDS_id].get_seq()[phase_slice_offset:]).translate() + '\n')
                        cds_out.close()
                    lengths_dict[cluster_index] = exon_lengths[:]
                    
def main():
    if sys.argv[1] == 'output_fastas':
        clusters = pickle.load(open('_clusters.pkl','rb'))
        prot_seq_dict = pickle.load(open('_psd.pkl','rb'))
        if '[' in sys.argv[4]:
            ref_genomes = eval(sys.argv[4])
            annotations = eval(sys.argv[5])
        else:
            ref_genomes = sys.argv[4]
            annotatios = sys.argv[5]
        output_fastas(clusters,prot_seq_dict,sys.argv[2],sys.argv[3],ref_genomes, annotations)

if __name__ == "__main__":
    main()

    
