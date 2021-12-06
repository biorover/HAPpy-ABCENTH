#!/usr/bin/env python
#Toolbox for HAP.py

#import genome_fork as genome
import subprocess
import time
import re

def genewise2aughints(genewiselines, src="G"):
    """Finds gff lines in a list of genewise lines and converts them to augustus hint format"""
    outlines = []
    for line in genewiselines:
        if line.count('\t') > 7:
            fields = line.split('\t')
            if fields[6] == '-':
                fields[3],fields[4] = fields[4],fields[3]
            if fields[2] == 'cds':
                outlines.append('\t'.join(fields[:2] + ['exon'] + fields[3:8] + ['src=' + src]))
            elif fields[2] == 'intron':
                if fields[6] == '+':
                    dss = fields[3]
                    ass = fields[4]
                elif fields[6] == '-':
                    dss = fields[4]
                    ass = fields[3]
                outlines.append('\t'.join(fields[:2] + ['dss',dss,dss] + fields[5:8] + ['src=' + src]))
                outlines.append('\t'.join(fields[:2] + ['ass',ass,ass] + fields[5:8] + ['src=' + src]))
    return(outlines)
    

def parse_genewise(genewise_file, use_gff_adjust = True, name_prefix = ""):
    """Pulls out gff lines and fasta sequences from genewise output"""
    gff_lines = []
    fasta_lines = []
    fasta_part = False
    for line in open(genewise_file):
        if line.count('\t') > 7 and "\tcds\t" in line:
            fields = line.split('\t')
            fields[8] = name_prefix + fields[8]
            fields[2] = "CDS"
            if use_gff_adjust:
                offset = int(fields[0].split('_')[1].split('-')[0])
                new_locus = fields[0].split('_')[0]
                gff_lines.append(gff_adjust(["\t".join(fields)],offset = offset,new_locus = new_locus)[0])
            else:
                gff_lines.append("\t".join(fields))
        elif line[0] == '>' and not "Results for " in line:
            fasta_lines.append(">" + name_prefix + line[1:])
            fasta_part = True
        elif line[0] == '/':
            fasta_part = False
        elif fasta_part:
            fasta_lines.append(line)
    return (gff_lines, fasta_lines)

def parse_augustus(augustus_file,name_prefix = None, use_gff_adjust = True):
    gtf_lines = []
    fasta_lines = []
    fasta_add = False
    for line in open(augustus_file):
        if line.count('\t') > 7 and "\tCDS\t" in line:
            fields = line.split('\t')
            working_name = fields[8].split('"')[-2]
            if name_prefix:
                working_name = name_prefix + "_" + working_name
                fields[8] = 'gene_id ' + working_name + ';transcript_id ' + working_name + '-RA\n'
            if use_gff_adjust:
                offset = int(fields[0].split('_')[1].split('-')[0])
                new_locus = fields[0].split('_')[0]
                gtf_lines.append(gff_adjust(["\t".join(fields)],offset,new_locus)[0])
            else:
                gtf_lines.append("\t".join(fields))
        elif "# protein sequence = [" in line:
            if not ']' in line:
                fasta_add = True
            fasta_lines.append('>' + working_name + '\n' + line.split('[')[1].replace(']',''))
        elif ']' in line and fasta_add:
            fasta_add = False
            fasta_lines.append(line[2:].replace(']',''))
        elif fasta_add:
            fasta_lines.append(line[2:])
    return (gtf_lines, fasta_lines)


def gff_adjust(gff_list,offset,new_locus = None):
    """Adds offset to each coordinate in gff_list"""
    adjusted_lines = []
    for line in range(len(gff_list)):
        fields = gff_list[line].split('\t')
        start,stop = fields[3],fields[4]
        fields[3] = str(min([int(start) + offset,int(stop) + offset]))
        fields[4] = str(max([int(start) + offset,int(stop) + offset]))
        if new_locus:
            fields[0] = new_locus  
        adjusted_lines.append('\t'.join(fields))
    return(adjusted_lines)

