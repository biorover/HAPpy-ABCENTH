#!/usr/bin/env python

import sys
import os
import shutil
import argparse
from HAPpy import genome_fork as genome
from HAPpy import toolbox_for_HAP
from HAPpy import CandidateLociBuilder
from HAPpy import AnnotatorRunner
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

def happy_logo():
    return """            _______  _______                    _
  |\     /|(  ___  )(  ____ )                  ( )
  | )   ( || (   ) || (    )|                  | |
  | (___) || (___) || (____)|                  | |
  |  ___  ||  ___  ||  _____) _______          | |
  | (   ) || (   ) || (      (  ____ )|\     /|(_)
  | )   ( || )   ( || )  _   | (    )|( \   / ) _
  |/     \||/     \||/  (_)  | (____)| \ (_) / (_)
                             |  _____)  \   /
                             | (         ) (
                             | )         | |
                             |/          \_/
\n"""

####
#
# accessory functions for multithreading
#
####

def run_genewisedb(genewisedb_path,hmmconvert_path,hmm3_file,seq_file,start,end,out_dir):
    tmp_fasta = tempfile.NamedTemporaryFile()
    tmp_hmm = tempfile.NamedTemporaryFile()
    tmp_out = tempfile.TemporaryFile()
    subprocess.call(shlex.split(hmmconvert_path + ' -2 ' + hmm3_file),stdout = tmp_hmm,stderr = open('/dev/null','w'))
    tmp_hmm.flush()
    seq_file_file = open(seq_file).read().split('\n')
    seqname = seq_file_file[0][1:]
    tmp_fasta.write(bytes('>' + seqname + '\n','UTF-8'))
    tmp_fasta.write(bytes("".join(seq_file_file[1:])[start:end],'UTF-8'))
    hmm_name = hmm3_file.split('/')[-1].replace('.hmm','')
    tmp_fasta.flush()
    genewise_result = subprocess.check_output(shlex.split(genewisedb_path +
                                ' -cut 10 -sum -gff -hmmer ' + tmp_hmm.name + ' ' + tmp_fasta.name),
                                stderr = open('/dev/null')).decode('utf-8')
    hit_name = hmm_name + '_hitOn_' + seqname + '_' + str(start) + '-' + str(end)
    if "GeneWise\tmatch\t" in genewise_result:
        out_file = open(out_dir + '/' + hit_name + '.genewise.gff','w')
        hit_file = open(out_dir + '/' + hit_name + '.genewise.hits','w')
        next_lines_are_summaries = False
        for line in genewise_result.split('\n'):
            fields = line.split('\t')
            if len(fields) > 8:
                if fields[2] == 'match':
                    score = fields[5]
                elif fields[2] in ['cds','intron']:
                    if fields[2] == 'cds':
                        fields[2] = 'CDS'
                    fields[5] = score
                    qstart = min([int(fields[3]),int(fields[4])])
                    qend = max([int(fields[3]),int(fields[4])])
                    fields[3] = str(qstart + start)
                    fields[4] = str(qend + start)
                    hit_suffix = fields[8].split('-')[-1]
                    fields[8] = 'gene_id ' + hit_name + '-' + hit_suffix + ';transcript_id ' + hit_name + '-' + hit_suffix + '-RA\n'
                    out_file.write("\t".join(fields))
            elif 'Bits' in line and 'Query' in line:
                next_lines_are_summaries = True
            elif "//" in line:
                next_lines_are_summaries = False
            elif next_lines_are_summaries:
                sfields = line.split()
                if float(sfields[0]) > 0:
                    if int(sfields[5]) > int(sfields[6]):
                        strand = "-"
                    else:
                        strand = "+"
                    hit_file.write("\t".join([sfields[1],sfields[4],".",".",".",strand,sfields[2],sfields[3],str(int(sfields[5]) + start),str(int(sfields[6]) + start),'.',sfields[0],'.','.']) + '\n')
        hit_file.close()
        out_file.close()

def run_exonerate(exonerate_path,query_fasta,seq_file,start,end,out_dir, exonerate_options,splice3,splice5):
    tmp_fasta = tempfile.NamedTemporaryFile()
    tmp_out = tempfile.TemporaryFile()
    seq_file_file = open(seq_file).read().split('\n')
    seqname = seq_file_file[0][1:]
    tmp_fasta.write(bytes('>' + seqname + '\n','UTF-8'))
    tmp_fasta.write(bytes("".join(seq_file_file[1:])[start:end],'UTF-8'))
    tmp_fasta.flush()
    query_name = query_fasta.split('/')[-1].replace('.fa','').replace('.consensus','')
    exonerate_result = subprocess.check_output(shlex.split(exonerate_path + ' ' + exonerate_options +
                                ' --showtargetgff -q ' + query_fasta + ' -t ' + tmp_fasta.name + ' --splice3 ' +
                                splice3 + ' --splice5 ' + splice5 )).decode('utf-8')
    hit_name = query_name + '_hitOn_' + seqname + '_' + str(start) + '-' + str(end)
    if "C4 Alignment:" in exonerate_result:
        out_file = open(out_dir + '/' + hit_name + '.exonerate.gff','w')
        hit_file = open(out_dir + '/' + hit_name + '.exonerate.hits','w')
        for line in exonerate_result.split('\n'):
            fields = line.split('\t')
            if len(fields) > 8:
                if fields[2] == 'gene':
                    score = fields[5]
                    hit_suffix = fields[8].split()[1]
                if fields[2] in ['cds','intron']:
                    if fields[2] == 'cds':
                        fields[2] = 'CDS'
                    fields[5] = score
                    qstart = min([int(fields[3]),int(fields[4])])
                    qend = max([int(fields[3]),int(fields[4])])
                    fields[3] = str(qstart + start)
                    fields[4] = str(qend + start)
                    fields[8] = 'gene_id ' + hit_name + '-' + hit_suffix + ';transcript_id ' + hit_name + '-' + hit_suffix + '-RA\n'
                    out_file.write("\t".join(fields))
            elif "vulgar:" in line:
                fields = line.split()
                hit_file.write("\t".join([query_name,fields[5],'.','.','.',fields[8],str(int(fields[2]) + 1),fields[3],str(int(fields[6]) + start),str(int(fields[7]) + start),'.',fields[9],'.','.']) + '\n')
        out_file.close()
        hit_file.close()

####
#
# main program functions
#
####

def set_program_paths(path_list = None):
    """sets program file paths using either default names or paths supplied with
    as a comma seperated list (e.g. via the "program_filepaths" argument)"""
    path_dict = {}
    path_dict['program_dir'] = "/".join(os.path.abspath(__file__).split('/')[:-1]) + '/'
    path_dict['hmmsearch'] = "hmmsearch"
    path_dict['hmmconvert'] = "hmmconvert"
    path_dict['hmmbuild'] = "hmmbuild"
    path_dict['mafft'] = "mafft"
    path_dict['thammerin'] = path_dict['program_dir'] + 'thammerin.py'
    path_dict['augustus'] = "augustus"
    path_dict['genewisedb'] = "genewisedb"
    path_dict['diamond'] = 'diamond'
    path_dict['hmmemit'] = 'hmmemit'
    path_dict['exonerate'] = 'exonerate'
    if path_list:
        for line in path_list.split(','):
            exec('path_dict["' + line.replace('=','"]="') + '"')
    return path_dict

def args_check(args):
    """checks run mode arguments to make sure they are compatible with how HAP runs"""
    if args.annotations:
        if not args.ref_genome:
             sys.exit('Argument error: A reference genome file must be provide under "--ref_genome" for each gtf given under "--annotations"')
        elif type(args.annotations) == list:
            if len(args.annotations) != len(args.ref_genome):
                sys.exit('Argument error: A reference genome file must be provide under "--ref_genome" for each gtf given under "--annotations"')
    if args.search_mode == "exons" and not args.annotations:
        sys.exit('Argument error: "--search_mode exons" can only be used with the input options "--annotations" + "--ref_genome"')
    if args.annotator == "ABCENTH":
        if args.search_mode != "exons" and not args.hmm_dir:
            sys.exit('Argument error: "--annotator ABCENTH" can only be used with "--search_mode exons" and the input options "--annotations" + "--ref_genome"')

def translate_genome(genome_fasta,out_file,min_orf_size,max_orf_size):
    target_nucdb = genome.Genome(genome_fasta, truncate_names = True)
    frame_fasta = open(out_file,'w')
    for seq_id in target_nucdb.genome_sequence:
        frameonef = genome.Sequence(target_nucdb.genome_sequence[seq_id])
        frameoner = genome.Sequence(target_nucdb.genome_sequence[seq_id]).reverse_compliment()
        frames = [frameonef,genome.Sequence(frameonef[1:]), genome.Sequence(frameonef[2:]),
                  frameoner, genome.Sequence(frameoner[1:]),genome.Sequence(frameoner[2:])]
        fasta_list = []
        for frame_num in (0,1,2,3,4,5):
            frame = frames[frame_num]
            if frame_num < 3:
                frame_offset = frame_num
            else:
                frame_offset = len(frameonef) - frame_num + 3
            orfs = frame.translate(trimX = False).split('*')
            last_orf_end = 0
            for orf in orfs:
                orf_start = last_orf_end
                last_orf_end = last_orf_end + 1 + len(orf)
                if len(orf) > min_orf_size and len(orf) < max_orf_size:
                    fasta_list.append('>' + seq_id + "_frameOffset-" + str(frame_offset) + "_orfStart-" + str(orf_start) +
                                      '\n' + orf)
        frame_fasta.write('\n'.join(fasta_list) + '\n')
    frame_fasta.close()

def annotations2fl(annotations,ref_genome,output_file):
    """outputs full length protein sequences from annotations for clustering"""
    fasta_out = open(output_file,'w')
    for i in range(len(annotations)):
        sp_genome = genome.Genome(ref_genome[i])
        sp_genome.read_gff(annotations[i])
        fasta_out.write(sp_genome.annotations.get_fasta('transcript',seq_type = 'protein') + "\n")
    fasta_out.close()

def build_clusters(protein_fasta,out_dir,threads,dthresh,path_dict,log_file):
    """Builds a pairwise-alignment UPGMA tree using MAFFT and breaks the tree into clusters of genes\
    within a specified p-distance of each other"""
    prot_seq_dict = {}
    for line in open(protein_fasta):
        if line[0] == ">":
            active_seq = line[1:-1]
            prot_seq_dict[active_seq] = ""
        else:
            prot_seq_dict[active_seq] += line[:-1]
    sys.stderr.write('building clusters for ' + str(len(prot_seq_dict)) + ' protiens in input fasta\n')
    subprocess.call(shlex.split(path_dict['mafft'] + ' --thread ' + str(threads) + ' --globalpair --treeout ' + protein_fasta),stdout = open(out_dir + '/flprots.mafft','w'), stderr = log_file)
    prottree = ete3.Tree(open(out_dir + '/flprots.fa.tree').read() + ';')
    clusters = []
    dont_append = False
    while prottree.get_distance(prottree.get_leaves()[0]) > dthresh and len(prottree.get_leaves()) > 1:
        children = prottree.get_children()
        child_one_dist = children[0].get_distance(children[0].get_leaves()[0])
        child_two_dist = children[1].get_distance(children[1].get_leaves()[0])
        while child_one_dist > dthresh and child_two_dist > dthresh:
            children = children[0].get_children()
            child_one_dist = children[0].get_distance(children[0].get_leaves()[0])
            child_two_dist = children[1].get_distance(children[1].get_leaves()[0])
        if child_one_dist < dthresh:
            clusters.append(children[0].get_leaf_names())
            try:
                prottree.prune(list(set(prottree.get_leaves()) - set(children[0].get_leaves())), preserve_branch_length = True)
            except:
                if len(prottree.get_leaves()) < 2:
                    dont_append = True
                    break
        if child_two_dist < dthresh:
            clusters.append(children[1].get_leaf_names())
            try:
                prottree.prune(list(set(prottree.get_leaves()) - set(children[1].get_leaves())), preserve_branch_length = True)
            except:
                if len(prottree.get_leaves()) < 2:
                    dont_append = True
                    break
    if not dont_append:
        clusters.append(prottree.get_leaf_names())
    return clusters,prot_seq_dict

def build_hmms(fasta_dir,out_dir,path_dict,mafft_options,threads,pre_aligned = False):
    """creates hidden markov models for every fasta in the "fasta_dir" using hmmbuild
    and writes these hmms to the "out_dir"."""
    fasta_files = os.listdir(fasta_dir)
    sys.stderr.write('aligning and building hmms for ' + str(len(fasta_files)) + ' fasta files\n')
    threads_list = []
    if not pre_aligned:
        for fasta in fasta_files:
            file_root = ".".join(fasta.split('.')[:-1])
            threads_list.append(
                threading.Thread(target = subprocess.call, args = [shlex.split(path_dict['mafft'] + " " + mafft_options + " " \
                                 + fasta_dir + "/" + fasta)],
                                 kwargs = {'stdout':open(fasta_dir + "/" + file_root + ".mafft",'w'),'stderr':open('/dev/null','w')})
            )
            threads_list[-1].start()
            while threading.active_count() >= threads:
                time.sleep(0.1)
        for thread in threads_list:
            thread.join()
    threads_list = []
    for fasta in fasta_files:
        file_root = ".".join(fasta.split('.')[:-1])
        if pre_aligned:
            align_file = fasta
        else:
            align_file = file_root + '.mafft'
        threads_list.append(
            threading.Thread(target = subprocess.call, args = [shlex.split(path_dict['hmmbuild'] + " --amino " + out_dir + "/" + file_root + ".hmm "\
                        + fasta_dir + '/' + align_file)],
                        kwargs = {'stderr':open('/dev/null','w'),'stdout':open('/dev/null','w')})
        )
        threads_list[-1].start()
        while threading.active_count() >= threads:
            time.sleep(0.1)
    for thread in threads_list:
        thread.join()

def filter_cluster_hmms_by_alignments(alignment_dir,hmm_dir,filter_options,filtered_hmm_dir):
    keep_cluster = []
    for alignment_file in os.list_dir(alignment_dir):
        if ".mafft" in alignment_file:
            pass

def run_thammerin(hmm_dir,target_file,thmmerin_dir,threads,path_dict,genome_orfs,out_dir,evalue):
    hmm_list = os.listdir(hmm_dir)
    threads_list = []
    sys.stderr.write('running ' + str(len(hmm_list)) + ' thammerin searches\n')
    decile_written = []
    for i in range(len(hmm_list)):
        hmm = hmm_list[i]
        decile = float( str(1.0 * i / len(hmm_list))[:3])
        if decile in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9] and not decile in decile_written:
            sys.stderr.write('    finished ' + str(int(100 * decile)) + '% of thammerin searches\n')
            decile_written.append(decile)
        file_root = ".".join(hmm.split('.')[:-1])
        if not genome_orfs:
            subprocess.call(shlex.split('thammerin --frames_out '
                            + out_dir + "/genomeFrames.fa -e " + str(evalue) +
                            " -p " + hmm_dir + '/' + hmm + ' -n ' + target_file +
                            ' -f ' + path_dict['hmmsearch']),
                            stdout = open(thmmerin_dir + '/' + file_root + '.tab','w'))
            genome_orfs = out_dir + '/genomeFrames.fa'
        else:
            threads_list.append(threading.Thread(target = subprocess.call,
                                args = [shlex.split('thammerin' +
                                " -e " + str(evalue) +
                                " -p " + hmm_dir + '/' + hmm + ' --frames_in ' + genome_orfs +
                                ' -f ' + path_dict['hmmsearch'])],
                                kwargs = {'stdout':open(thmmerin_dir + '/' + file_root + '.tab','w')}
            ))
            threads_list[-1].start()
        while threading.active_count() >= 2 + threads / 3:
            time.sleep(0.1)
    for thread in threads_list:
        thread.join()
    thammerin_hits = open(out_dir + '/thammerin_hits.tsv','w')
    for thammerin_file in os.listdir(thmmerin_dir):
        for line in open(thmmerin_dir + '/' + thammerin_file):
            thammerin_hits.write(line)
    thammerin_hits.close()

def emit_hmm_consensi(hmm_dir,out_dir,path_dict):
    hmm_files = os.listdir(hmm_dir)
    for hmm in hmm_files:
        out_file = open(out_dir + '/' + hmm.replace('.hmm','') + '.consensus.fa','w')
        emittext = subprocess.check_output(shlex.split(path_dict['hmmemit'] + ' -c ' + hmm_dir + '/' +
                                    hmm)).decode('utf-8')
        out_file.write(emittext.replace('-consensus',''))
        out_file.close()

def run_diamond(query_peps,diamonddb,path_dict,threads,out_file_name,diamond_options,evalue):
    diamond_hits = subprocess.check_output(shlex.split(path_dict['diamond'] + ' blastp ' + diamond_options + ' -f 6 sseqid qseqid pident length mismatch gapopen sstart send qstart qend evalue bitscore qlen slen --query '
            + query_peps + ' --db ' + diamonddb + ' --evalue ' +  str(evalue) + ' --threads ' + str(threads)),stderr = open('/dev/null','w') ).decode('utf-8')
    out_file = open(out_file_name,'w')
    for line in diamond_hits.split('\n'):
        fields = line.split('\t')

        if len(fields) > 10:
            tname = fields[1]
            frame_offset = int(tname.split('_frameOffset-')[1].split('_orfStart-')[0])
            orf_start = int(tname.split('_orfStart-')[1])
            if frame_offset < 3:
                strand = "+"
                start = orf_start * 3 + frame_offset + int(fields[8]) * 3 - 2
                stop = orf_start * 3 + frame_offset + int(fields[9]) * 3
            else:
                strand = "-"
                start = frame_offset - 3 * orf_start - 3 * int(fields[8]) + 3
                stop = frame_offset - 3 * orf_start - 3 * int(fields[9]) + 1
            out_file.write("\t".join([fields[0],tname.split('_frameOffset-')[0],fields[2],'.','.',strand,
                             fields[6],fields[7],str(start),str(stop),fields[10],fields[11],fields[12],fields[13]]) + '\n')
    out_file.close()

def hit_table2candidate_loci(hit_table,search_mode,max_loci_per_cluster,max_intron):
    candidate_loci_dict = {}
    locus_min_scores = {}
    for line in open(hit_table):
        fields = line.split('\t')
        if search_mode == 'exons':
            cluster = fields[0].split('exon')[0]
        else:
            cluster = fields[0]
        if not cluster in candidate_loci_dict:
            candidate_loci_dict[cluster] = {}
            locus_min_scores[cluster] = []
        # parses lines from table
        seqname = fields[1]
        qstart = int(fields[6])
        qend = int(fields[7])
        tcoords = [int(fields[8]),int(fields[9])]
        tstart = min(tcoords)
        tend = max(tcoords)
        score = float(fields[11])
        strand = fields[5]
        #
        # here we will add each overlapping hit to a list, then go back and modify them to merge/extend as necessary
        overlapped_hits = []
        new_locus = [tstart,tend,qstart,qend,strand,score]
        if not seqname in candidate_loci_dict[cluster]: #no hits yet for this cluster/sequence, can add new one
            candidate_loci_dict[cluster][seqname] = []
        else:
            left_append = True
            right_append = True
            for locus_num in range(len(candidate_loci_dict[cluster][seqname])):
                locus = candidate_loci_dict[cluster][seqname][locus_num]
                if strand == locus[4] and locus[0] - max_intron < tstart < tend < locus[1] + max_intron and (
                        (tstart < locus[0] and qstart < locus[2] and strand == "+" and left_append) or
                        (tend > locus[1] and qend > locus[3] and strand == "+" and right_append) or
                        (tstart < locus[0] and qend > locus[3] and strand == "-" and left_append) or
                        (tend > locus[1] and qstart < locus[2] and strand == "-" and right_append) or
                        locus[0] < tstart < tend < locus[1] ):
                    overlapped_hits.append(locus_num)
                    if tstart < locus[0]:
                        left_append = False
                    if tend > locus[1]:
                        right_append = False
                    if locus[0] < new_locus[0]:
                        new_locus[0] = locus[0]
                    if locus[1] > new_locus[1]:
                        new_locus[1] = locus[1]
                    if locus[2] < new_locus[2]:
                        new_locus[2] = locus[2]
                    if locus[3] > new_locus[3]:
                        new_locus[3] = locus[3]
                    new_locus[5] += locus[5]
        overlapped_hits.reverse() # done so that iterating over this to delete will delete highest numbers first
        for locus_num in overlapped_hits:
            del candidate_loci_dict[cluster][seqname][locus_num]
        candidate_loci_dict[cluster][seqname].append(new_locus)
        candidate_loci_dict[cluster][seqname].sort()
        locus_min_scores[cluster].append(new_locus[5])
    if max_loci_per_cluster < float('inf'):
        for cluster in locus_min_scores:
            score_copy = locus_min_scores[cluster][:]
            score_copy.sort()
            locus_min_scores[cluster] = score_copy[-1 * int(max_loci_per_cluster):]
    candidate_locus_list = []
    for cluster in candidate_loci_dict:
        for seqname in candidate_loci_dict[cluster]:
            for locus in candidate_loci_dict[cluster][seqname]:
                if locus[5] >= min(locus_min_scores[cluster]):
                    candidate_locus_list.append([cluster,seqname,str(locus[0]),str(locus[1]),locus[4],str(locus[2]),str(locus[3]),str(locus[5])])
    return candidate_locus_list

def filter_canidate_loci(candidate_loci):
    """sorts candidate loci by score (highest first) and then filters out (or trims) loci which overlap a higher scoring locus"""
    seq_trees = {}
    sorted_loci = []
    for fields in candidate_loci:
        score = float(fields[7])
        sorted_loci.append([score,fields])
    sorted_loci.sort()
    sorted_loci.reverse()
    for entry in sorted_loci:
        fields = entry[1]
        seq = fields[1]
        strand = fields[4]
        if not seq + strand in seq_trees:
            seq_trees[seq + strand] = IntervalTree()
        start = int(fields[2])
        end = int(fields[3])
        new_coords = [[start]]
        completely_overlapped = False
        for overlap in sorted(seq_trees[seq+strand][start:end]):
            if overlap[0] <= start and overlap[1] < end:
                new_coords[-1][0] = overlap[1] + 1
            elif overlap[1] < end:
                new_coords[-1].append(overlap[0])
                new_coords.append([overlap[1] + 1])
            elif overlap[1] >= end and overlap[0] > start:
                new_coords[-1].append(overlap[0])
            elif overlap[1] >= end and overlap[0] <= start:
                completely_overlapped = True
        if len(new_coords[-1]) < 2:
            new_coords[-1].append(end)
        if not completely_overlapped:
            for coords in new_coords:
                if coords[1] - coords[0] > 100:
                    newfields = fields[:]
                    newfields[2] = str(coords[0])
                    newfields[3] = str(coords[1])
                    seq_trees[seq + strand][coords[0]:coords[1]] = newfields
    filtered_candidate_loci = []
    for tree in seq_trees:
        for entry in sorted(seq_trees[tree]):
            filtered_candidate_loci.append(entry[2])
    return filtered_candidate_loci

def annotate_with_genewise(genome_file,buffer,candidate_loci_file,path_dict,hmm_dir,out_dir,threads):
    sys.stderr.write('annotating with genewise\n')
    genewisedb_path = path_dict['genewisedb']
    hmmconvert_path = path_dict['hmmconvert']
    out_file = open('/dev/null')
    for line in open(genome_file):
        if line[0] == ">":
            out_file.close()
            current_seq = line[1:].split()[0]
            out_file = open(out_dir + '/' + current_seq + '.fasta','w')
        out_file.write(line)
    out_file.close()
    threads_list = []
    for line in open(candidate_loci_file):
        fields = line.split()
        hmm3_file = hmm_dir + '/' + fields[0] + '.hmm'
        start = max([1,int(fields[2]) - buffer])
        end = int(fields[3]) + buffer
        seq_file = out_dir + '/' + fields[1] + '.fasta'
        threads_list.append(threading.Thread(target = run_genewisedb,
            args = [genewisedb_path,hmmconvert_path,hmm3_file,seq_file,start,end,out_dir]))
        threads_list[-1].start()
        while threading.active_count() > threads:
            time.sleep(0.1)
    for thread in threads_list:
        thread.join()
    out_file = open(out_dir + '/all_genewise_predictions.gff','w')
    for genewise_file in os.listdir(out_dir):
        if 'genewise.gff' in genewise_file:
            for line in open(out_dir + '/' + genewise_file):
                out_file.write(line)
    target_genome = genome.Genome(genome_file)
    target_genome.read_gff(out_dir + '/all_genewise_predictions.gff',features_to_ignore = ['intron'])
    pep_out = open(out_dir + '/all_genewise_predictions.pep','w')
    pep_out.write(target_genome.annotations.get_fasta('gene',seq_type='protein'))
    pep_out.close()
    out_file.close()

def annotate_with_exonerate(genome_file,buffer,candidate_loci_file,path_dict,fasta_dir,out_dir,threads,splice3,splice5,exonerate_options):
    sys.stderr.write('annotating with exonerate\n')
    exonerate_path = path_dict['exonerate']
    out_file = open('/dev/null')
    for line in open(genome_file):
        if line[0] == ">":
            out_file.close()
            current_seq = line[1:].split()[0]
            out_file = open(out_dir + '/' + current_seq + '.fasta','w')
        out_file.write(line)
    out_file.close()
    threads_list = []
    for line in open(candidate_loci_file):
        fields = line.split('\t')
        seq_file = out_dir + '/' + fields[1] + '.fasta'
        query_fasta = fasta_dir + '/' + fields[0] + '.consensus.fa'
        start = max([int(fields[2]) - buffer,0])
        end = int(fields[3]) + buffer
        threads_list.append(threading.Thread(target = run_exonerate,
            args = [exonerate_path,query_fasta,seq_file,start,end,out_dir, exonerate_options,splice3,splice5] ))
        threads_list[-1].start()
        while threading.active_count() > threads:
            time.sleep(0.1)
    out_file = open(out_dir + '/all_exonerate_predictions.gff','w')
    for thread in threads_list:
        thread.join()
        thread.join()
    for exonerate_file in os.listdir(out_dir):
        if 'exonerate.gff' in exonerate_file:
            for line in open(out_dir + '/' + exonerate_file):
                out_file.write(line)
    out_file.close()
    target_genome = genome.Genome(genome_file)
    target_genome.read_gff(out_dir + '/all_exonerate_predictions.gff',features_to_ignore = ['intron'])
    pep_out = open(out_dir + '/all_exonerate_predictions.pep','w')
    pep_out.write(target_genome.annotations.get_fasta('gene',seq_type='protein'))
    pep_out.close()

def merge_hits(out_file,*args):
    seqstranddict = {}
    out = open(out_file,'w')
    for hit_file in args:
        for line in open(hit_file):
            fields = line.split('\t')
            if len(fields) > 10:
                seqstrand = fields[1] + fields[5]
                if not seqstrand in seqstranddict:
                    seqstranddict[seqstrand] = IntervalTree()
                tcoords = [int(fields[8]),int(fields[9])]
                qcoords = [int(fields[6]),int(fields[7])]
                tcoords.sort()
                score = float(fields[11])
                target = fields[0]
                dellist = []
                for overlap in sorted(seqstranddict[seqstrand][tcoords[0]:tcoords[1]]):
                    if overlap[2][0] == target:
                        tcoords = [min(tcoords + overlap[2][1]),max(tcoords + overlap[2][1])]
                        qcoords = [min(qcoords + overlap[2][3]),max(qcoords + overlap[2][3])]
                        score = max([score,overlap[2][2]])
                        dellist.append(overlap)
                for overlap in dellist:
                    seqstranddict[seqstrand].discard(overlap)
                if seqstrand[-1] == '-':
                    write_tcoords = [tcoords[1],tcoords[0]]
                else:
                    write_tcoords = tcoords
                fields[8] = str(write_tcoords[0])
                fields[9] = str(write_tcoords[1])
                fields[11] = str(score)
                seqstranddict[seqstrand][tcoords[0]:tcoords[1]] = [target,tcoords,score,qcoords,"\t".join(fields)]
    for seqstrand in seqstranddict:
        for interval in sorted(seqstranddict[seqstrand]):
            out.write(interval[2][4])
    out.close()

def annotate_with_augustus(genome_file,augustus_species,user_hints,profile_dir,hmm_hints,out_dir,candidate_loci_file,buffer,config_file,genewise,exonerate,path_dict,threads,filter_by_overlap,trim_neighboring_loci):
    sys.stderr.write('annotating candidate loci with augustus\n')
    log_file = open(out_dir + '/augustus_cmds.log','w')
    err_log_file = open(out_dir + '/augustus_errors.log','w')
    out_file = open('/dev/null')
    for line in open(genome_file):
        if line[0] == ">":
            out_file.close()
            current_seq = line[1:].split()[0]
            out_file = open(out_dir + '/' + current_seq + '.fasta','w')
        out_file.write(line)
    out_file.close()
    if user_hints:
        hints_files = {}
        for line in open(user_hints):
            seq = line.split()[0]
            if not seq in hints_files:
                hints_files[seq] = open(seq + '.hints.gff','w')
            hints_files[seq].write(line)
        for hint_file in hints_files:
            hints_files[hint_file].close()
    if hmm_hints:
        hints_files = {}
        for line in open(hmm_hints):
            fields = line.split('\t')
            seq = fields[1]
            cluster = fields[0]
            hint_start = min([int(fields[9]),int(fields[8])]) + 10
            hint_end = max([int(fields[9]),int(fields[8])]) - 10
            hint_file_name = seq + '_' + cluster
            if hint_end - hint_start > 5:
                hint_out = open(out_dir + '/' + hint_file_name + '.hints.gff','a')
                hint_out.write("\t".join([seq,'thammerin','CDSpart',str(hint_start),str(hint_end),fields[11],fields[5],'.',
                                        'grp=' + cluster + '-searchHit;pri=4;src=P\n']))
            hint_out.close()
    if genewise:
        hints_files = {}
        for line in open(genewise):
            fields = line.split('\t')
            seq = fields[0]
            hint_start = fields[3]
            hint_end = fields[4]
            score = fields[5]
            strand = fields[6]
            target = fields[8].split(';')[0].split('_hitOn_')[0].split()[1]
            feature = fields[2]
            hint_file_name = seq + '_' + target
            hint_out = open(out_dir + '/' + hint_file_name + '.hints.gff','a')
            hint_out.write("\t".join([seq,'GeneWise',feature,hint_start,hint_end,score,strand,'.',
                                    'grp=' + target + '-genewise;pri=2;src=P\n']))
            hint_out.close()
    if exonerate:
        for line in open(exonerate):
            fields = line.split('\t')
            seq = fields[0]
            hint_start = fields[3]
            hint_end = fields[4]
            score = fields[5]
            strand = fields[6]
            target = fields[8].split(';')[0].split('_hitOn_')[0].split()[1]
            feature = fields[2]
            hint_file_name = seq + '_' + target
            if feature == 'intron' and int(hint_end) - int(hint_start) < 100: #add this because augustus throws out hint groups with short introns
                continue
            hint_out = open(out_dir + '/' + hint_file_name + '.hints.gff','a')
            hint_out.write("\t".join([seq,'exonerate',feature,hint_start,hint_end,score,strand,'.',
                                    'grp=' + target + '-exonerate;pri=2;src=P\n']))
            hint_out.close()
    aug_cmd_root = path_dict['augustus'] + ' --species=' + augustus_species + ' --alternatives-from-evidence=false --alternatives-from-sampling=false --stopCodonExcludedFromCDS=false'
    threads_list = []
    loci_dict = {}
    for line in open(candidate_loci_file):
        fields = line.split('\t')
        if not fields[1] + fields[4] in loci_dict:
            loci_dict[fields[1] + fields[4]] = []
        loci_dict[fields[1] + fields[4]].append([int(fields[2]),fields])
    last_locus = ['',0]
    for seq in loci_dict: # sets annotation start and end points, adjusting them so that the annotation stards and ends (locus starts and ends +- buffer) don't overlap with previous locus starts and ends (unless ends are less than 100bp from each other, in which case they can overlap by up to 100bp)
        loci_dict[seq].sort()
        for locus_index in range(len(loci_dict[seq])):
            locus = loci_dict[seq][locus_index]
            fields = locus[1]
            start = max([1,int(fields[2]) - buffer])
            end = int(fields[3]) + buffer
            if last_locus[0] == seq and start <= last_locus[1] and trim_neighboring_loci: #edit the following lines if I decide to make annotation regions completely non-overlapping (currently they are just filtered to not overlap the previous/next hit range)
                start = min([last_locus[1] + 1, int(fields[2]) - 100])
                loci_dict[seq][locus_index - 1][3] = max([int(fields[2]) - 1, loci_dict[seq][locus_index - 1][3] - buffer + 100])
            loci_dict[seq][locus_index] += [start,end]
            last_locus = [seq,int(fields[3])]
    for seqstrand in loci_dict:
        for locus in loci_dict[seqstrand]:
            fields = locus[1]
            start = locus[2]
            end = locus[3]
            aug_cmd = aug_cmd_root + ' --predictionStart=' + str(start) + \
                ' --predictionEnd=' + str(int(fields[3]) + buffer)
            if seqstrand[-1] == '+':
                aug_cmd += ' --strand=forward'
            elif seqstrand[-1] == '-':
                aug_cmd += ' --strand=backward'
            seq = fields[1]
            query = fields[0]
            if seq + '.hints.gff' in os.listdir(out_dir) or seq + '_' + query + '.hints.gff' in os.listdir(out_dir):
                if seq + '.hints.gff' in os.listdir(out_dir) and seq + '_' + query + '.hints.gff' in os.listdir(out_dir):
                    cl_hints_file = tempfile.NamedTemporaryFile('.gff')
                    for line in open(out_dir + '/' + seq + '.hints.gff','b'):
                        cl_hints_file.write(line)
                    for line in  open(out_dir + '/' + seq + '_' + query + '.hints.gff','b'):
                        cl_hints_file.write(line)
                    cl_hints_file.flush()
                elif seq + '.hints.gff' in os.listdir(out_dir):
                    cl_hints_file = out_dir + '/' + seq + '.hints.gff'
                else:
                    cl_hints_file = out_dir + "/" + seq + '_' + query + '.hints.gff'
                aug_cmd += ' --hintsfile=' + cl_hints_file + ' --extrinsicCfgFile=' + config_file
            if profile_dir:
                aug_cmd += ' --proteinprofile=' + profile_dir + '/' + query + '.prfl'
            aug_cmd += " " + out_dir + '/' + seq + '.fasta'
            out_file_name = query + '_seqID_' + seq + '_coords_' + fields[2] + "-" + fields[3] + '.aug.gff'
            log_file.write(aug_cmd + '\n')
            threads_list.append(threading.Thread(target = subprocess.call,
                args = [shlex.split(aug_cmd)], kwargs = {'stdout':open(out_dir + '/' + out_file_name,'a'),
                                                        'stderr':err_log_file}
            ))
            threads_list[-1].start()
            while threading.active_count() > threads:
                time.sleep(0.1)
    for thread in threads_list:
        thread.join()
    master_augustus_file = open(out_dir + '/all_augustus_predictions.gff','w')
    if filter_by_overlap:
        overlap_tree_dict = {}
        for line in open(filter_by_overlap):
            fields = line.split('\t')
            start = min([int(fields[9]),int(fields[8])])
            end = max([int(fields[9]),int(fields[8])])
            cluster = fields[0]
            seq = fields[1]
            if not cluster in overlap_tree_dict:
                overlap_tree_dict[cluster] = {}
            if not seq in overlap_tree_dict[cluster]:
                overlap_tree_dict[cluster][seq] = IntervalTree()
            overlap_tree_dict[cluster][seq][start:end] = True
    for out_file in os.listdir(out_dir):
        if ".aug.gff" in out_file:
            gene_id_root = out_file.replace('aug.gff','')
            cluster = gene_id_root.split('_seqID_')[0]
            seq = gene_id_root.split('_seqID_')[1].split('_coords_')[0]
            can_loc_start = int(gene_id_root.split('_seqID_')[1].split('_coords_')[1].split('-')[0])
            can_loc_end =  int(gene_id_root.split('_seqID_')[1].split('_coords_')[1].split('-')[1].split('.')[0])
            if filter_by_overlap:
                outlinedict = {}
                goodgenes = set()
            for line in open(out_dir + '/' + out_file):
                if "AUGUSTUS\tCDS\t" in line:
                    if filter_by_overlap:
                        overlap_tree = IntervalTree(overlap_tree_dict[cluster][seq][can_loc_start:can_loc_end])
                        fields = line.split('\t')
                        if len(overlap_tree[int(fields[3]):int(fields[4])]) > 0:
                            goodgenes.add(fields[8])
                        if not fields[8] in outlinedict:
                            outlinedict[fields[8]] = []
                        outlinedict[fields[8]].append(line.replace(' "g',' "' + gene_id_root))
                    else:
                        master_augustus_file.write(line.replace(' "g',' "' + gene_id_root))
            if filter_by_overlap:
                for goodgene in goodgenes:
                    master_augustus_file.write("".join(outlinedict[goodgene]))
    master_augustus_file.close()
    target_genome = genome.Genome(genome_file)
    target_genome.read_gff(out_dir + '/all_augustus_predictions.gff')
    pep_out = open(out_dir + '/all_augustus_predictions.pep','w')
    pep_out.write(target_genome.annotations.get_fasta('gene',seq_type='protein'))
    pep_out.close()
    log_file.close()
    err_log_file.close()

def main():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                 usage = "\npython HAP.py [optional arguments] --genome <genome.fa> " +
                                 " > output.gtf\npython HAP.py [optional arguments] --genome <genome.fa> ".join([
                                 "--annotations <ann.gtf> [ann2.gtf ...] --ref_genome <ref.fa> [ref2.fa ...]",
                                 "--protein_seqs <proteins.fa>", "--hmm <proteins.hmm>", "--fasta_dir <fasta_directory/>",
                                 "--alignment_dir <alignment_directory/>","--hmm_dir <hmm_directory/> [--augustus_profile_dir <prfl_directory>]"])
                                 + ' > output.gtf',
                                 description= happy_logo() + "Pipeline for annotating genes in a genome using homologous sequences \
    from a related species. \n\nProgram dependencies: \
    python, mafft, hmmer suite v3, ete3 python library, and genewise. \n\nHAP.py can \
    be run with a single set of related genes, multiple predefined clusters of \
    related genes, or it can build clusters from large highly divergent gene families. \
    HAP.py can take as input a gtf and genome from one or more related species, \
    an unaligned fasta file of query protein sequences, an HMM built from query \
    protein sequences using HMMER v3, or directories containing files corresponding to pre-defined \
    clusters of related genes- either as unaligned fasta files, protein sequence alignments, or HMMs.")

    manditory_args = parser.add_argument_group(title = 'mandatory arguments')
    query_args = parser.add_argument_group(title = "query arguments (at least one required)")
    addl_args = parser.add_argument_group(title = "additional HAP.py arguments")
    run_args = parser.add_argument_group(title = "client program arguments")
    manditory_args.add_argument('--genome',dest='target_genome', help = 'genome to be annotated')
    query_args.add_argument('--annotations', nargs = "*", default = None, help = 'One or more sets of annotations (in gtf format; requires "--ref_genome")')
    query_args.add_argument('--ref_genome', nargs = "*", default = None, help = 'One or more genomes of closely related species\
                        (you should provide one genome for each gtf given for --annotations and ensure the order is the same; genomes should \
                        be in fasta format)')
    query_args.add_argument('--thammerin_results', dest = 'thammerin_results', default = None, help = 'pre-computed thammerin results')
    query_args.add_argument('--protein_seqs', default = None,
                        help = 'homologous protein sequences (can replace -a + -r for the "full_length" runmode)')
    query_args.add_argument('--hmm', default = None, help = "Protein HMM (built with HMMER v3 hmmbuild)")
    query_args.add_argument('--fasta_dir', default = None, help = 'Directory of UNALIGNED fasta files for predefined clusters')
    query_args.add_argument('--alignment_dir', default = None, help = 'Directory of ALIGNED fasta files for predefined clusters')
    query_args.add_argument('--hmm_dir', default = None, help = "Directory of HMM files for predefined clusters")
    query_args.add_argument('--hit_table', default = None, help = 'precomputed hit table')
    run_args.add_argument('--program_filepaths', default = None,
                        help = 'optional comma seperated list of file paths for programs not in PATH variable, formated as "programName1=/path/to/program1,programName2=/path/to/program2"')
    addl_args.add_argument('--min_orf_size', default = 10, type = int,
                        help = 'minimum size for orfs to be search by hmmsearch')
    addl_args.add_argument('--cutoff', default = 1.0, type = float,
                        help = 'Distance cutoff for seperating proteins into clusters. Accepts values from zero to one, default = 1.0 (no breaking into clusters)')
    addl_args.add_argument('--threads', default = 1, type = int,
                        help = "number of threads to be used with processes that support multithreading (mafft and thammerin)")
    addl_args.add_argument('--buffer', dest = 'buffer', default = 5000, type = int,
                        help = 'buffer on either side of loci identified to feed into gene predictor (default = 5000)')
    addl_args.add_argument('--evalue', default = 0.01, type = float, help = 'Evalue cutoff for thammerin')
    addl_args.add_argument('--genome_orfs', default = None, help = 'ORFs file from previous thammerin run (saves about five minutes for insect-sized genomes)')
    addl_args.add_argument('--search_mode', default = 'fl', help = 'Search with full length sequences ("fl") or individual exons ("exons"). \
                        Default = "fl". "exons" requires input in the form of --annotations + --ref_genome, and this option is required for \
                        the "--annotator ABCENTH" option.')
    addl_args.add_argument('--annotator', default = "genewise", help = 'Program to use for building final annotations. Currently the options \
                        are "genewise" (default) and "ABCENTH". I plan to add support for hint-guided AUGUSTUS at some point.')
    run_args.add_argument('--ABCENTH_full_pseudoexon_search', default = True, help = 'Speficies whether ABCENTH should use genewise to \
                          extend pseudogenized exon hits (default = True). Specifying "False" can substantially speed up ABCENTH on highly fragmented genome \
                          at the expense of recovering slightly less sequence for pseudogenes (genes without frameshifts/internal stop codons are unaffected)')
    run_args.add_argument('--augustus_profile_dir', default = None, help = 'directory of augustus profiles corresponding to entries in "--hmm_dir" \
                        (for optional use with "--hmm_dir <dir>" and "--annotator augustus")')
    run_args.add_argument('--augustus_species', default = 'fly', help = 'species name for augustus parameters (default = "fly")')
    addl_args.add_argument('--output_dir', default = 'HAPpy_results',help = 'folder to write results to (default = "HAPpy_results")')
    addl_args.add_argument('--overwrite', default = False, type = bool, help = 'overwrite specified output dir (default = False)')
    run_args.add_argument('--cluster_mafft_options', default = '--globalpair --maxiterate 1000', help = 'options for aligning sequences \
                        within clusters using mafft (default = "--globalpair --maxiterate 1000")')
    addl_args.add_argument('--max_loci_per_cluster', default = float('inf') , type = float,
                            help = 'maximum number of candidate loci to annotate per cluster (default = inf)')
    addl_args.add_argument('--max_intron_length', default = 20000,type = int, help = 'maximum length of potential introns (used to decide if hits from\
                        the same cluster on the same sequence are a single locus or two seperate loci; default = 200000)')
    run_args.add_argument('--augustus_extrinsic', default = None, help = 'augustus extrinsic config file')
    run_args.add_argument('--augustus_hints', default = None, help = 'gff file of hints to pass to augustus for annotation')
    run_args.add_argument('--augustus_use_hmm_hints', default = True, help = 'Pass thmmerin hits to augustus for annotation')
    addl_args.add_argument('--initial_search_tool',default = 'thammerin',help = 'Search tool for initial identification of candidate loci. \
                            Accepts "thammerin" (default) or "diamond".')
    run_args.add_argument('--diamond_options',default = "--very-sensitive", help = 'Additional arguments to pass to the diamond aligner (if using diamond instead of thammerin). \
                        default = "--very-sensitive"')
    run_args.add_argument('--exonerate_options',default = '--model protein2genome --percent 10',
                            help = 'Additional arguments to pass to exonerate protein2genome (if using exonerate for annotation). \
                            default = "--model protein2genome --percent 10". For exhaustive exonerate annotation, we recomend "--model protein2genome:bestfit -E -S no --percent 20"')
    addl_args.add_argument('--cluster_filters',nargs = "*",default = None, help = 'prefilter for clusters before taking them forward in analysis (default = None). \
                        provide space-seperated list of "param=value" Keys. Valid params: "min_num_species" (requires species prefix folowed by underscore in gene name), \
                        "max_alignment_gap_percent" (for total alignment), "max_gaps_for_min_species" (min_num_species must also be set, at least min_num_species \
                        must have no more than this percentage of gaps in alignment)')
    addl_args.add_argument('--extend_loci_with_alignments',default = "True",type = eval, help = 'if exonerate or genewise is run and augustus is to be run, use alignments from exonerate \
                        and/or genewise to potentially extend candidate loci before running augustus (default = True)')
    addl_args.add_argument('--trim_neighboring_loci', default = "True", type = eval, help = 'Prevents augustus candidate gene prediction regions (candidate loci += buffer) from overlapping neighboring candidate loci (NB: buffers \
                        can still overlap each other when True, core candidate loci still won\'t overlap if False) (default = True)')

    args = parser.parse_args()
    sys.stderr.write(happy_logo())
    sys.stderr.write('called as "' + " ".join(sys.argv) + '"\n')
    path_dict = set_program_paths(args.program_filepaths)
    args_check(args)
    if os.path.isdir(args.output_dir):
        if args.overwrite:
            subprocess.call(['rm','-r',args.output_dir])
        else:
            sys.exit('output directory exists, and "overwrite" argument not set to True')
    os.makedirs(args.output_dir + '/clusters')
    os.makedirs(args.output_dir + '/hmms')
    os.makedirs(args.output_dir + '/thammerin_results')
    run_log = open(args.output_dir + '/happy_run_log.out','w')
    run_log.write('HAPp called as "' + " ".join(sys.argv) + '"\n')
    if args.annotations:
        annotations2fl(args.annotations,args.ref_genome,args.output_dir + '/flprots.fa')
    if args.protein_seqs:
        shutil.copyfile(args.protein_seqs,args.output_dir + '/flprots.fa')
    if args.protein_seqs or args.annotations:
        if args.cutoff < 1.0:
            clusters,prot_seq_dict = build_clusters(args.output_dir + '/flprots.fa',
                                                    args.output_dir,
                                                    args.threads,args.cutoff, path_dict,run_log)
            pickle.dump(clusters,open('_clusters.pkl','wb'))
            pickle.dump(prot_seq_dict,open('_psd.pkl','wb'))
            subprocess.call(['python', '-m', 'HAPpy.HAPpy-hiMem-client.py','output_fastas',args.search_mode,
                          args.output_dir + '/clusters',str(args.ref_genome), str(args.annotations)])
            del clusters
            del prot_seq_dict
        else:
            shutil.copyfile(args.output_dir + '/flprots.fa',args.output_dir + '/clusters/flprots.fa')
    if args.fasta_dir:
        fasta_dir = args.fasta_dir
    else:
        fasta_dir = args.output_dir + '/clusters'
    if args.fasta_dir or args.protein_seqs or args.annotations:
        build_hmms(fasta_dir,args.output_dir + '/hmms',path_dict,args.cluster_mafft_options,args.threads,False)
    elif args.alignment_dir:
        build_hmms(args.alignment_dir,args.output_dir + '/hmms',path_dict,args.cluster_mafft_options,args.threads,True)
    elif args.hmm:
        shutil.copy(args.hmm,args.output_dir + '/hmms/')
    if args.hmm_dir:
        hmm_dir = args.hmm_dir
    elif args.hmm or args.fasta_dir or args.protein_seqs or args.annotations:
        hmm_dir = args.output_dir + '/hmms'
    else:
        hmm_dir = None
    if hmm_dir and (args.initial_search_tool == 'diamond' or 'exonerate' in args.annotator):
        os.makedirs(args.output_dir + '/consensus_fastas')
        emit_hmm_consensi(hmm_dir,args.output_dir + '/consensus_fastas',path_dict)
    if (args.hmm_dir or args.hmm or args.fasta_dir or args.protein_seqs or args.annotations) and not args.hit_table:
        if not args.genome_orfs:
            translate_genome(args.target_genome,args.output_dir + '/genomeOrfs.fasta',args.min_orf_size,10000)
            genome_orfs = args.output_dir + '/genomeOrfs.fasta'
        else:
            genome_orfs = args.genome_orfs
        if args.initial_search_tool == 'thammerin':
            run_thammerin(hmm_dir,args.target_genome,args.output_dir + '/thammerin_results',
                      args.threads,path_dict,genome_orfs,args.output_dir,args.evalue)
            hit_table = args.output_dir + '/thammerin_hits.tsv'
        elif args.initial_search_tool == 'diamond':
            consensus_peps = open(args.output_dir + '/all_query_consensi.pep','w')
            for fasta_file in os.listdir(args.output_dir + '/consensus_fastas'):
                for line in open(args.output_dir + '/consensus_fastas/' + fasta_file):
                    consensus_peps.write(line)
            consensus_peps.close()
            subprocess.call(shlex.split(path_dict['diamond'] + ' makedb --in ' + args.output_dir +
                '/all_query_consensi.pep --db ' + args.output_dir + '/all_query_consensi' + ' --threads ' + str(args.threads) ),
                stderr = open('/dev/null'),stdout = open('/dev/null'))
            run_diamond(genome_orfs,args.output_dir + '/all_query_consensi',path_dict,args.threads,
                        args.output_dir + '/diamond_hits.tsv',args.diamond_options,args.evalue)
            hit_table = args.output_dir + '/diamond_hits.tsv'
    elif args.hit_table:
        hit_table = args.hit_table
    subprocess.call(shlex.split('sort -k2,2 -k9,9n ' + hit_table),stdout = open(args.output_dir + '/sorted_hit_table.tsv','w'))
    candidate_loci = hit_table2candidate_loci(args.output_dir + '/sorted_hit_table.tsv',args.search_mode,args.max_loci_per_cluster,args.max_intron_length)
    unfiltered_can_locs = open(args.output_dir + '/unfiltered_candidate_loci.tab','w') #adding this really just to debug
    for line in candidate_loci:
        unfiltered_can_locs.write("\t".join(line) + '\n')
    unfiltered_can_locs.close()
    filtered_candidate_loci = filter_canidate_loci(candidate_loci)
    candidate_locus_file = open(args.output_dir + '/candidate_loci.tab','w')
    for line in filtered_candidate_loci:
        candidate_locus_file.write("\t".join(line) + '\n')
    candidate_locus_file.close()
    if "ABCENTH" in args.annotator:
        subprocess.call(shlex.split('python -m ABCENTH.ParseHAPpyTable --table ' + args.output_dir + '/sorted_hit_table.tsv'),
            stdout = open(args.output_dir + '/abcenth_hit_table.tsv','w'))
        subprocess.call(shlex.split('python -m ABCENTH.HitTabFilter --table ' + args.output_dir + '/abcenth_hit_table.tsv'),
            stdout = open(args.output_dir + '/abcenth_filtered_hit_table.tsv','w'))
        subprocess.call(shlex.split('python -m ABCENTH.ParseHAPpyTable --hmm_dir '  + hmm_dir),
            stdout = open(args.output_dir + '/abcenth_cluster_info.tsv','w'))
        subprocess.call(shlex.split('ABCENTH --genome ' + args.target_genome + ' --table ' +
            args.output_dir + '/abcenth_filtered_hit_table.tsv --query_exon_info_table ' + args.output_dir + 
            '/abcenth_cluster_info.tsv --full_pseudoexon_search ' + str(args.ABCENTH_full_pseudoexon_search)), stdout = open(args.output_dir + '/ABCENTH.gtf','w'))
        
    if 'genewise' in args.annotator:
        os.environ['WISECONFIGDIR'] = path_dict['program_dir'] + '/' + 'cfgFiles'
        os.makedirs(args.output_dir + '/genewise')
        annotate_with_genewise(args.target_genome,args.buffer,args.output_dir + '/candidate_loci.tab',
            path_dict,hmm_dir,args.output_dir + '/genewise',args.threads)
    if 'exonerate' in args.annotator:
        os.makedirs(args.output_dir + '/exonerate/')
        splice3,splice5 = path_dict['program_dir'] + '/' + 'cfgFiles/spliceGT.txt',path_dict['program_dir'] + '/' + 'cfgFiles/spliceAG.txt'
        annotate_with_exonerate(args.target_genome,args.buffer,args.output_dir + '/candidate_loci.tab',
            path_dict,args.output_dir + '/consensus_fastas/',args.output_dir + '/exonerate/',args.threads,splice3,splice5,args.exonerate_options)
    if args.extend_loci_with_alignments and 'augustus' in args.annotator:
        merge_hits_args = [args.output_dir + '/merged_hits.tsv',args.output_dir + '/sorted_hit_table.tsv']
        if 'exonerate' in args.annotator:
            for exonerate_file in os.listdir(args.output_dir + '/exonerate/'):
                if "exonerate.hits" in exonerate_file:
                    merge_hits_args.append(args.output_dir + '/exonerate/' + exonerate_file)
        if 'genewise' in args.annotator:
            for genewise_file in os.listdir(args.output_dir + '/genewise/'):
                if "genewise.hits" in genewise_file:
                    merge_hits_args.append(args.output_dir + '/genewise/' + genewise_file)
        merge_hits(*merge_hits_args)
        candidate_loci = hit_table2candidate_loci(args.output_dir + '/merged_hits.tsv',args.search_mode,args.max_loci_per_cluster,args.max_intron_length)
        filtered_candidate_loci = filter_canidate_loci(candidate_loci)
        candidate_locus_file = open(args.output_dir + '/alignment_extened_candidate_loci.tab','w')
        for line in filtered_candidate_loci:
            candidate_locus_file.write("\t".join(line) + '\n')
        candidate_locus_file.close()
        augustus_candidate_loci = args.output_dir + '/alignment_extened_candidate_loci.tab'
    else:
        augustus_candidate_loci = args.output_dir + '/candidate_loci.tab'
    if "augustus" in args.annotator:
        os.makedirs(args.output_dir + '/augustus')
        if args.augustus_use_hmm_hints:
            hmm_hints = hit_table
        else:
            hmm_hints = None
        if args.augustus_extrinsic:
            extrinsic_config_file = args.augustus_extrinsic
        else:
            extrinsic_config_file = path_dict['program_dir'] + 'cfgFiles/augustus_extrinsic.MP.cfg'
        if 'genewise' in args.annotator:
            genewise = args.output_dir + '/genewise/all_genewise_predictions.gff'
        else:
            genewise = False
        if 'exonerate' in args.annotator:
            exonerate = args.output_dir + '/exonerate/all_exonerate_predictions.gff'
        else:
            exonerate = False
        annotate_with_augustus(args.target_genome,args.augustus_species,args.augustus_hints,
            args.augustus_profile_dir,hmm_hints,args.output_dir + '/augustus',
            augustus_candidate_loci,args.buffer,
            extrinsic_config_file,genewise,exonerate,path_dict,args.threads,hit_table,args.trim_neighboring_loci)

if __name__ == "__main__":
    main()


####
#
# features to add:
# -genewise result filtering (by what?)
# -make limit on minimum candidate locus length
# -function to generate augustus profiles
# -evaluation of cluster matches (number, length, "goodness", etc.)
# -automated cluster selection by conservation (select relatively gapless, conserved clusters for continued proccessing to add in annotation efforts)
# -automatic retraining and re-running of augustus
# -install script
#
# Stretch goals:
# -function to check for and report frameshifts using genewise/exonerate matches
# -gene-location aware candidate locus parsing
# -phase-aware genewise analysis
# -rewrite to avoid using thammerin
# -put on conda
#
# Finished:
# -output peptides from exonerate, genewise and augustus analyses
# -evidence modeler makes predictions substantially worse, so this project was abandoned
# -automatic extension of candidate loci based on genewise/exonerate hits
#
####
