#!/usr/bin/env python
import argparse
import os
parser = argparse.ArgumentParser(description='splits query name output by HAP.py and builds table required for ABCENTH')
parser.add_argument('--table',default = None, help = 'table output by HAP.py')
parser.add_argument('--hmm_dir',default = None, help = "director with all cluster hmms")
args = parser.parse_args()

if args.table:
    for line in open(args.table):
        fields = line.replace('\n','').replace('\r','').split('\t')
        cluster = fields[0].split('exon')[0]
        exon_number = fields[0].split('exon')[1].split('of')[0]
        number_of_exons = fields[0].split('of')[1].split('phases')[0]
        start_phase = fields[0].split('phases')[1].split('and')[0]
        end_phase = fields[0].split('and')[1].split('.')[0]
        aa_len = fields[12]
        print('\t'.join([cluster] + fields[1:12] + [start_phase,end_phase,aa_len,exon_number,number_of_exons]))

elif args.hmm_dir:
    for hmm_file in os.listdir(args.hmm_dir):
        if hmm_file[-4:] == ".hmm" and not "fullLenForHMM" in hmm_file:
            cluster = hmm_file.split('exon')[0]
            exon_number = hmm_file.split('exon')[1].split('of')[0]
            number_of_exons = hmm_file.split('of')[1].split('phases')[0]
            start_phase = hmm_file.split('phases')[1].split('and')[0]
            end_phase = hmm_file.split('and')[1].split('.')[0]
            aa_len = open(args.hmm_dir + "/" + hmm_file).read().split('\n')[2].split()[1].replace('\r','')
            print('\t'.join([cluster,exon_number,number_of_exons,start_phase,end_phase,aa_len,os.path.abspath(args.hmm_dir) + '/' + hmm_file]))
