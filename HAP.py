#!/usr/bin/env python

####
#
# features to add:
# -gene-location aware candidate locus parsing
# -phase-aware genewise analysis
#
####

import sys
import os
import shutil
import argparse
import genome_fork as genome
import toolbox_for_HAP
import CandidateLociBuilder
import AnnotatorRunner
import subprocess
import ete3
import time
import copy
import threading
import shlex
import tempfile



parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                 usage = "\npython HAP.py [optional arguments] --genome <genome.fa> " +
                                 " > output.gtf\npython HAP.py [optional arguments] --genome <genome.fa> ".join([
                                 "--annotations <ann.gtf> [ann2.gtf ...] --ref_genome <ref.fa> [ref2.fa ...]",
                                 "--protein_seqs <proteins.fa>", "--hmm <proteins.hmm>", "--fasta_dir <fasta_directory/>",
                                 "--alignment_dir <alignment_directory/>","--hmm_dir <hmm_directory/> [--augustus_profile_dir <prfl_directory>]"])
                                 + ' > output.gtf',
                                 description= """            _______  _______                    _
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
  \n""" + "Pipeline for annotating genes in a genome using homologous sequences \
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
run_args.add_argument('--program_filepaths', default = None,
                    help = 'optional comma seperated list of file paths for programs not in PATH variable, formated as "programName1=/path/to/program1,programName2=/path/to/program2"')
addl_args.add_argument('--min_orf_size', default = 10, type = int,
                    help = 'minimum size for orfs to be search by hmmsearch')
addl_args.add_argument('--cutoff', default = 1.0, type = float,
                    help = 'Distance cutoff for seperating proteins into clusters. Accepts values from zero to one, default = 1.0 (no breaking into clusters)')
addl_args.add_argument('--threads', default = 1, type = int,
                    help = "number of threads to be used with processes that support multithreading (mafft and thammerin)")
query_args.add_argument('--thammerin_results', dest = 'thammerin_results', default = None, help = 'pre-computed thammerin results')
addl_args.add_argument('--buffer', dest = 'buffer', default = 5000, type = int,
                    help = 'buffer on either side of loci identified to feed into gene predictor (default = 5000)')
addl_args.add_argument('--evalue', default = 0.01, type = float, help = 'Evalue cutoff for thammerin')
addl_args.add_argument('--genome_orfs', default = None, help = 'ORFs file from previous thammerin run (saves about five minutes for insect-sized genomes)')
addl_args.add_argument('--search_mode', default = 'fl', help = 'Search with full length sequences ("fl") or individual exons ("exons"). \
                    Default = "fl". "exons" requires input in the form of --annotations + --ref_genome, and this option is required for \
                    the "--annotator ABCENTH" option.')
addl_args.add_argument('--annotator', default = "genewise", help = 'Program to use for building final annotations. Currently the options \
                    are "genewise" (default) and "ABCENTH". I plan to add support for hint-guided AUGUSTUS at some point.')
query_args.add_argument('--protein_seqs', default = None,
                    help = 'homologous protein sequences (can replace -a + -r for the "full_length" runmode)')
query_args.add_argument('--hmm', default = None, help = "Protein HMM (built with HMMER v3 hmmbuild)")
query_args.add_argument('--fasta_dir', default = None, help = 'Directory of UNALIGNED fasta files for predefined clusters')
query_args.add_argument('--alignment_dir', default = None, help = 'Directory of ALIGNED fasta files for predefined clusters')
query_args.add_argument('--hmm_dir', default = None, help = "Directory of HMM files for predefined clusters")
run_args.add_argument('--augustus_profile_dir', default = None, help = 'directory of augustus profiles corresponding to entries in "--hmm_dir" \
                    (for optional use with "--hmm_dir <dir>" and "--annotator augustus")')
run_args.add_argument('--augustus_species', default = 'fly', help = 'species name for augustus parameters (default = "fly")')
addl_args.add_argument('--output_dir', default = 'HAPpy_results',help = 'folder to write results to (default = "HAPpy_results")')
query_args.add_argument('--hit_table', default = None, help = 'precomputed hit table')
addl_args.add_argument('--overwrite', default = False, type = bool, help = 'overwrite specified output dir (default = False)')
run_args.add_argument('--cluster_mafft_options', default = '--globalpair --maxiterate 1000', help = 'options for aligning sequences \
                    within clusters using mafft (default = "--globalpair --maxiterate 1000")')
addl_args.add_argument('--max_loci_per_cluster', default = float('inf') , type = float,
                        help = 'maximum number of candidate loci to annotate per cluster (default = inf)')
addl_args.add_argument('--max_intron_length', default = 20000,help = 'maximum length of potential introns (used to decide if hits from\
                       the same cluster on the same sequence are a single locus or two seperate loci; default = 200000)')
run_args.add_argument('--augustus_extrinsic', default = None, help = 'augustus extrinsic config file')
run_args.add_argument('--augustus_hints', default = None, help = 'gff file of hints to pass to augustus for annotation')
run_args.add_argument('--augustus_use_hmm_hints', default = True, help = 'Pass thmmerin hits to augustus for annotation')

args = parser.parse_args()

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
    seqname = seq_file_file[0][1:-1]
    tmp_fasta.write(bytes('>' + seqname + '\n','UTF-8'))
    tmp_fasta.write(bytes("".join(seq_file_file[1:])[start:end],'UTF-8'))
    #char = seq_file_file.read(1)
    #counter = 1
    hmm_name = hmm3_file.split('/')[-1].replace('.hmm','')
    #while char:
    #    char = seq_file_file.read(1)
    #    if counter >= start:
    #        tmp_fasta.write(bytes(char,'UTF-8'))
    #    if counter > end:
    #        break
    #    if char != '/n':
    #        counter += 1
    tmp_fasta.flush()
    genewise_result = subprocess.check_output(shlex.split(genewisedb_path +
                                ' -cut 10 -gff -hmmer ' + tmp_hmm.name + ' ' + tmp_fasta.name),
                                stderr = open('/dev/null')).decode('utf-8')
    hit_name = hmm_name + '_hitOn_' + seqname + '_' + str(start) + '-' + str(end)
    #print(genewise_result)
    if "GeneWise\tmatch\t" in genewise_result:
        out_file = open(out_dir + '/' + hit_name + '.genewise.gff','w')
        for line in genewise_result.split('\n'):
            fields = line.split('\t')
            if len(fields) > 8:
                if fields[2] == 'match':
                    score = fields[5]
                if fields[2] in ['cds','intron']:
                    if fields[2] == 'cds':
                        fields[2] = 'CDS'
                    fields[5] = score
                    qstart = min([int(fields[3]),int(fields[4])])
                    qend = max([int(fields[3]),int(fields[4])])
                    fields[3] = str(qstart + start)
                    fields[4] = str(qend + start)
                    hit_suffix = fields[8].split('-')[-1][:-1]
                    fields[8] = 'gene_id ' + hit_name + '-' + hit_suffix + ';transcript_id ' + hit_name + '-' + hit_suffix + '-RA\n'
                    out_file.write("\t".join(fields))
        out_file.close()

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
    if path_list:
        for line in path_list.split(','):
            exec("path_dict[" + line.replace('=',']='))
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
        if args.search_mode != "exons":
            sys.exit('Argument error: "--annotator ABCENTH" can only be used with "--search_mode exons" and the input options "--annotations" + "--ref_genome"')

def annotations2fl(annotations,ref_genome,output_file):
    """outputs full length protein sequences from annotations for clustering"""
    fasta_out = open(output_file,'w')
    for i in range(len(annotations)):
        sp_genome = genome.Genome(ref_genome[i])
        sp_genome.read_gff(annotations[i])
        fasta_out.write(sp_genome.annotations.get_fasta('transcript',seq_type = 'protein'))
    fasta_out.close()

def build_clusters(protein_fasta,out_dir,threads,dthresh,path_dict):
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
    subprocess.call(shlex.split(path_dict['mafft'] + ' --thread ' + str(threads) + ' --globalpair --treeout ' + out_dir + '/flprots.fa'),stdout = open(out_dir + '/flprots.mafft','w'), stderr = open('/dev/null'))
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

def output_fastas(clusters, prot_seq_dict, search_mode, fasta_dir,ref_genome_list = None,annotations = None):
    """Uses clusters dictionary built by "build_clusters" to break input sequences into fastas\
    for individual clusters. If search_mode = "exons", builds fastas for each exon."""
    sys.stderr.write('writing fastas for ' + str(len(clusters)) + ' clusters\n')
    cluster_dict = {}
    lengths_dict = {}
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
                ref_genomes = []
                for i in range(len(ref_genome_list)):
                    ref_genomes.append(genome.Genome(ref_genome_list[i]))
                    ref_genomes[-1].read_gff(annotations[i])
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
                            sys.stder.write("uh oh, transcripts with varying number of exons in: cluster " + str(cluster_index) + ", specifically gene " + transcript_id)
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
            subprocess.call(shlex.split(path_dict['thammerin'] + ' --frames_out '
                            + out_dir + "/genomeFrames.fa -e " + str(evalue) +
                            " -p " + hmm_dir + '/' + hmm + ' -n ' + target_file +
                            ' -f ' + path_dict['hmmsearch']),
                            stdout = open(thmmerin_dir + '/' + file_root + '.tab','w'))
            genome_orfs = out_dir + '/genomeFrames.fa'
        else:
            threads_list.append(threading.Thread(target = subprocess.call,
                                args = [shlex.split(path_dict['thammerin'] +
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

def annotate_with_augustus(genome_file,augustus_species,user_hints,profile_dir,hmm_hints,out_dir,candidate_loci_file,buffer,config_file,genewise,path_dict,threads):
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
            hint_start = min([int(fields[9]),int(fields[8])]) + 10
            hint_end = max([int(fields[9]),int(fields[8])]) - 10
            if hint_end - hint_start > 5:
                if not seq in hints_files:
                    hints_files[seq] = open(out_dir + '/' + seq + '.hints.gff','a')
                hints_files[seq].write("\t".join([seq,'thammerin','CDSpart',str(hint_start),str(hint_end),fields[11],'.','.',
                                        'grp=' + fields[0] + ';pri=4;src=P\n']))
        for hint_file in hints_files:
            hints_files[hint_file].close()
    if genewise:
        hints_files = {}
        for line in open(genewise):
            fields = line.split('\t')
            seq = fields[0]
            hint_start = fields[3]
            hint_end = fields[4]
            score = fields[6]
            target = fields[8].split(';')[0].split('_hitOn_')[0].split()[1]
            feature = fields[2]
            if not seq in hints_files:
                hints_files[seq] = open(out_dir + '/' + seq + '.hints.gff','a')
            hints_files[seq].write("\t".join([seq,'GeneWise',feature,hint_start,hint_end,score,'.','.',
                                    'grp=' + target + ';pri=2;src=P\n']))
        for hint_file in hints_files:
            hints_files[hint_file].close()
    aug_cmd_root = path_dict['augustus'] + ' --protein=T --species=' + augustus_species
    threads_list = []
    for line in open(candidate_loci_file):
        fields = line.split('\t')
        aug_cmd = aug_cmd_root + ' --predictionStart=' + str(max([0,int(fields[2]) - buffer])) + \
            ' --predictionEnd=' + str(int(fields[3]) + buffer)
        seq = fields[1]
        query = fields[0]
        if seq + '.hints.gff' in os.listdir(out_dir):
            aug_cmd += ' --hintsfile=' + out_dir + '/' + seq + '.hints.gff' + ' --extrinsicCfgFile=' + config_file
        if profile_dir:
            aug_cmd += ' --proteinprofile=' + profile_dir + '/' + query + '.prfl'
        aug_cmd += " " + out_dir + '/' + seq + '.fasta'
        out_file_name = query + '_' + seq + '_' + fields[2] + "-" + fields[3] + '.aug.gff'
        log_file.write(aug_cmd + '\n')
        threads_list.append(threading.Thread(target = subprocess.call,
            args = [shlex.split(aug_cmd)], kwargs = {'stdout':open(out_dir + '/' + out_file_name,'a'),
                                                    'stderr':err_log_file}
        ))
        threads_list[-1].start()
        while threading.active_count() >= threads:
            time.sleep(0.1)
    for thread in threads_list:
        thread.join()
    master_augustus_file = open(out_dir + '/all_augustus_predictions.gff','w')
    for out_file in os.listdir(out_dir):
        if ".aug.gff" in out_file:
            for line in open(out_dir + '/' + out_file):
                master_augustus_file.write(line)
    master_augustus_file.close()
    log_file.close()
    err_log_file.close()

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
        while threading.active_count() >= threads:
            time.sleep(0.1)
    for thread in threads_list:
        thread.join()
    out_file = open(out_dir + '/all_genewise_predictions.gff','w')
    for genewise_file in os.listdir(out_dir):
        if 'genewise.gff' in genewise_file:
            for line in open(out_dir + '/' + genewise_file):
                out_file.write(line)
    out_file.close()


def main(args):
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
    sys.stderr.write("""            _______  _______                    _
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
\n""")
    if args.annotations:
        annotations2fl(args.annotations,args.ref_genome,args.output_dir + '/flprots.fa')
    if args.protein_seqs:
        shutil.copyfile(args.protein_seqs,args.output_dir + '/flprots.fa')
    if args.protein_seqs or args.annotations:
        if args.cutoff < 1.0:
            clusters,prot_seq_dict = build_clusters(args.output_dir + '/flprots.fa',
                                                    args.output_dir,
                                                    args.threads,args.cutoff, path_dict)
            output_fastas(clusters,prot_seq_dict,args.search_mode,
                          args.output_dir + '/clusters', args.ref_genome, args.annotations)
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
    if (args.hmm_dir or args.hmm or args.fasta_dir or args.protein_seqs or args.annotations) and not args.hit_table:
        run_thammerin(hmm_dir,args.target_genome,args.output_dir + '/thammerin_results',
                  args.threads,path_dict,args.genome_orfs,args.output_dir,args.evalue)
    if args.hit_table:
        hit_table = args.hit_table
    else:
        thammerin_hits = open(args.output_dir + '/thammerin_hits.tsv','w')
        for thammerin_file in os.listdir(args.output_dir + '/thammerin_results'):
            for line in open(args.output_dir + '/thammerin_results/' + thammerin_file):
                thammerin_hits.write(line)
        thammerin_hits.close()
        hit_table = args.output_dir + '/thammerin_hits.tsv'
    candidate_loci = hit_table2candidate_loci(hit_table,args.search_mode,args.max_loci_per_cluster,args.max_intron_length)
    candidate_locus_file = open(args.output_dir + '/candidate_loci.tab','w')
    for line in candidate_loci:
        candidate_locus_file.write("\t".join(line) + '\n')
    candidate_locus_file.close()
    if 'genewise' in args.annotator:
        os.environ['WISECONFIGDIR'] = path_dict['program_dir'] + '/' + 'happy_wisecfg'
        os.makedirs(args.output_dir + '/genewise')
        annotate_with_genewise(args.target_genome,args.buffer,args.output_dir + '/candidate_loci.tab',
            path_dict,hmm_dir,args.output_dir + '/genewise',args.threads)
    if "augustus" in args.annotator:
        os.makedirs(args.output_dir + '/augustus')
        if args.augustus_use_hmm_hints:
            hmm_hints = hit_table
        else:
            hmm_hints = None
        if args.augustus_extrinsic:
            extrinsic_config_file = args.augustus_extrinsic
        else:
            extrinsic_config_file = path_dict['program_dir'] + '/augustus_extrinsic.MP.cfg'
        if 'genewise' in args.annotator:
            genewise = args.output_dir + '/genewise/all_genewise_predictions.gff' #Fix after adding genewise function
        else:
            genewise = False
        annotate_with_augustus(args.target_genome,args.augustus_species,args.augustus_hints,
            args.augustus_profile_dir,hmm_hints,args.output_dir + '/augustus',
            args.output_dir + '/candidate_loci.tab',args.buffer,
            extrinsic_config_file,genewise,path_dict,args.threads)

if __name__ == "__main__":
    main(args)
