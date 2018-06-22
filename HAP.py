#!/usr/bin/env python

import sys
import os
import argparse
import genome_fork as genome
import toolbox_for_HAP
import CandidateLociBuilder
import AnnotatorRunner
import subprocess
import ete2
import time
import copy
import threading

parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                 usage = "\npython HAP.py [optional arguments] --genome <genome.fa> " +
                                 " > output.gtf\npython HAP.py [optional arguments] --genome <genome.fa> ".join([
                                 "--annotations <ann.gtf> [ann2.gtf ...] --ref_genome <ref.fa> [ref2.fa ...]",
                                 "--protein_seqs <proteins.fa>", "--hmm <proteins.hmm>", "--fasta_dir <fasta_directory/>",
                                 "--alignment_dir <alignment_directory/>","--hmm_dir <hmm_directory/> [--augustus_profile_dir <prfl_directory>]"])
                                 + ' > output.gtf',
                                 description= "Pipeline for annotating genes in a genome using homologous sequences \
from a related species. \n\nProgram dependencies: \
python, mafft, hmmer suite v3, ete2 python library, and genewise. \n\nHAP.py can \
be run with a single set of related genes, multiple predefined clusters of \
related genes, or it can build clusters from large highly divergent gene families. \
HAP.py can take as input a gtf and genome from one or more related species, \
an unaligned fasta file of query protein sequences, an HMM built from query \
protein sequences using HMMER v3, or directories containing files corresponding to pre-defined \
clusters of related genes- either as unaligned fasta files, protein sequence alignments, or HMMs.")

parser.add_argument('--genome',dest='target_genome', help = 'genome to be annotated')
parser.add_argument('--annotations', nargs = "*", default = None, help = 'One or more sets of annotations (in gtf format)')
parser.add_argument('--ref_genome', nargs = "*", default = None, help = 'One or more genomes of closely related species\
                    (you should provide one genome for each gtf given for --annotations and ensure the order is the same; genomes should \
                    be in fasta format)')
parser.add_argument('--program_filepaths', default = None, 
                    help = 'optional list of file paths for programs not in PATH variable, formated as "programName=/path/to/programName"')
parser.add_argument('--min_orf_size', default = 10, type = int,
                    help = 'minimum size for orfs to be search by hmmsearch')
parser.add_argument('--cutoff',dest = 'distance_cutoff', default = 1.0, type = float,
                    help = 'Distance cutoff for seperating proteins into clusters. Accepts values from zero to one, default = 1.0 (no breaking into clusters)')
parser.add_argument('--threads', default = 1, type = int,
                    help = "number of threads to be used with processes that support multithreading (mafft and thammerin)")
parser.add_argument('--thammerin_results', dest = 'thammerin_results', default = None, help = 'pre-computed thammerin results')
parser.add_argument('--buffer', dest = 'buffer', default = 500, type = int,
                    help = 'buffer on either side of loci identified to feed into gene predictor')
parser.add_argument('--evalue', default = 0.05, type = float, help = 'Evalue cutoff for thammerin')
parser.add_argument('--genome_orfs', default = None, help = 'ORFs file from previous thammerin run (saves about five minutes for insect-sized genomes)')
parser.add_argument('--search_mode', default = 'fl', help = 'Search with full length sequences ("fl") or individual exons ("exons"). \
                    Default = "fl". "exons" requires input in the form of --annotations + --ref_genome, and this option is required for \
                    the "--annotator ABCENTH" option.')
parser.add_argument('--annotator', default = "genewise", help = 'Program to use for building final annotations. Currently the options \
                    are "genewise" (default) and "ABCENTH". I plan to add support for hint-guided AUGUSTUS at some point.')
parser.add_argument('--protein_seqs', default = None,
                    help = 'homologous protein sequences (can replace -a + -r for the "full_length" runmode)')
parser.add_argument('--hmm', default = None, help = "Protein HMM (built with HMMER v3 hmmbuild)")
parser.add_argument('--fasta_dir', default = None, help = 'Directory of UNALIGNED fasta files for predefined clusters')
parser.add_argument('--alignment_dir', default = None, help = 'Directory of ALIGNED fasta files for predefined clusters')
parser.add_argument('--hmm_dir', default = None, help = "Directory of HMM files for predefined clusters")
parser.add_argument('--augustus_profile_dir', default = None, help = 'directory of augustus profiles corresponding to entries in "--hmm_dir" \
                    (for optional use with "--hmm_dir <dir>" and "--annotator augustus")')
parser.add_argument('--augustus_species', default = 'fly', help = 'species name for augustus parameters (default = "fly")')

args = parser.parse_args()

#sets up default program paths, overwritten by any program paths passes with -f or --program_filepaths
program_dir = "/".join(os.path.abspath(__file__).split('/')[:-1]) + '/'
hmmsearch = "hmmsearch"
hmmconvert = "hmmconvert"
hmmbuild = "hmmbuild"
mafft = "mafft"
thammerin = program_dir + 'thammerin.py'
augustus = "augustus"
genewisedb = "genewisedb"

if args.program_filepaths:
    for line in args.program_filepaths:
        exec(line.replace('\n','').replace('\r',''))


#checks that arguments make sense
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


#sets up some basic variables
dthresh = args.distance_cutoff
target_genome = os.path.abspath(args.target_genome)
call_directory = os.getcwd()
if args.augustus_profile_dir:
    augustus_profile_dir = os.path.abspath(args.augustus_profile_dir)
elif not args.hmm_dir:
    augustus_profile_dir = 'tmpOrAnPipePrfls'
else:
    augustus_profile_dir = None

if args.genome_orfs:
    thammerin_orfs = os.path.abspath(args.genome_orfs)
else:
    thammerin_orfs = 'tmpOrAnPipeFrames.fa'
if args.hmm_dir:
    hmm_dir = os.path.abspath(args.hmm_dir) + '/'
if args.fasta_dir:
    fasta_dir = os.path.abspath(args.fasta_dir) + '/'
if args.alignment_dir:
    alignment_dir = os.path.abspath(args.alignment_dir) + '/'
if args.protein_seqs:
    protein_seqs = os.path.abspath(args.protein_seqs)
    prot_seq_dict = {}
if args.hmm:
    protein_hmm = os.path.abspath(args.hmm)

#Reads in annotations, sets annotation genome path
if args.annotations:
    ref_genomes = []
    genome_paths = args.ref_genome
    geneset_paths = args.annotations
    for genome_index in range(len(genome_paths)):
        ref_genomes.append(genome.Genome(genome_paths[genome_index]))
        ref_genomes[genome_index].read_gff(geneset_paths[genome_index])


os.mkdir('tmpOrAnPipeDir')
os.chdir('tmpOrAnPipeDir')
os.mkdir('tmpOrAnPipeHMMs')
os.mkdir('tmpOrAnPipePrfls')

#Define functions which run individual steps of HAP.py program
def build_clusters():
    """Builds a pairwise-alignment UPGMA tree using MAFFT and breaks the tree into clusters of genes\
    within a specified p-distance of each other"""
    flprots = open('tmpOrAnPipe.flpeps.fa','w')
    if args.protein_seqs:
        protein_sequences = open(call_directory + '/' + args.protein_seqs).read()
        flprots.write(protein_sequences)
        for seq in protein_sequences.split('>')[1:]:
            seq_lines = seq.replace('\r','').split('\n')
            prot_seq_dict[seq_lines[0]] = "".join(seq_lines[1:])
    elif args.annotations:
        for ref_genome in ref_genomes:
            flprots.write(ref_genome.annotations.get_fasta('gene',seq_type = 'protein')+'\n')
    flprots.close()
    subprocess.call(mafft + ' --thread ' + str(args.threads) + ' --globalpair --treeout tmpOrAnPipe.flpeps.fa > tmpeOrAnPipe.tmp', shell = True)
    prottree = ete2.Tree(open('tmpOrAnPipe.flpeps.fa.tree').read() + ';')
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
    return clusters

def output_fastas(search_mode):
    """Uses clusters dictionary built by "build_clusters" to break input sequences into fastas\
    for individual clusters. If search_mode = "exons", builds fastas for each exon."""
    cluster_dict = {}
    lengths_dict = {}
    for cluster_index in range(len(clusters)):
        cluster = clusters[cluster_index]
        exon_file_root = "tmpOrAnPipeHMMs/cluster" + str(cluster_index)
        cluster_dict[cluster_index] = []
        exon_number = None
        exon_phases = None
        exon_lengths = None
        for seqname in cluster:   
            transcript_id = "_".join(seqname.split("_")[1:])
            if args.protein_seqs:
                seq_out = open(exon_file_root + '.fa','a')
                seq_out.write(">" + transcript_id + "\n" + prot_seq_dict[transcript_id] + '\n')
                seq_out.close()
            else:
                for ref_genome in ref_genomes:
                    try:
                        cluster_dict[cluster_index].append(ref_genome.annotations.transcript[transcript_id])
                    except KeyError:
                        continue
                    seq_out = open(exon_file_root + 'fullLenForHMM.fa','a')
                    seq_out.write(ref_genome.annotations.transcript[transcript_id].get_fasta(seq_type='protein') + '\n')
                    seq_out.close()
                    if search_mode == "exons":
                        if not exon_number:
                            exon_number = len(ref_genome.annotations.transcript[transcript_id].child_list)
                        else:
                            if len(ref_genome.annotations.transcript[transcript_id].child_list) != exon_number:
                                print "uh oh, transcripts with varying number of exons in: cluster " + str(cluster_index) + ", specifically gene " + transcript_id
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

def thammerin_func(file_name, hmm_dir):
    if not args.hmm_dir and not args.hmm:
        hmm_suffix = ".hmm"
        if open(query_pep_dir + file_name).read().count('>') > 1 and not args.alignment_dir:
            print "Aligning " + file_name
            subprocess.call("mafft --globalpair --maxiterate 1000 " + query_pep_dir + file_name +
                            " > tmpOrAnPipeHMMs/"+ file_name + ".mafftGinsi.fa 2> tmp.tmp",
                            shell = True)
        else:
            print "Copying " + file_name
            subprocess.call('cp ' + query_pep_dir + file_name + ' tmpOrAnPipeHMMs/' + file_name +
                            ".mafftGinsi.fa", shell = True)
        gene_cluster = file_name.split('exon')[0]
        if "augustus" in args.annotator and augustus_profile_dir:
            subprocess.call('msa2prfl.pl  tmpOrAnPipeHMMs/' + file_name + '.mafftGinsi.fa > tmpOrAnPipePrfls/'
                            + file_name + '.prfl', shell = True)
        try:
            subprocess.call(hmmbuild + " --amino tmpOrAnPipeHMMs/" + file_name + ".hmm tmpOrAnPipeHMMs/"
                            + file_name + ".mafftGinsi.fa > tmp.tmp", shell = True)
        except subprocess.CalledProcessError:
            print "Hm, problem with hmmbuild for file " + file_name
    else:
        hmm_suffix = ""
    query_len = int(subprocess.check_output('grep "LENG " ' + hmm_dir + file_name + hmm_suffix,
                                            shell = True).split()[1])
    subprocess.call("sed -e  '/NAME /s/$/_hmmFile_" + file_name + hmm_suffix + ".lenAdded.hmm" + "_len_" + str(query_len) + "/' " + hmm_dir + file_name +
                    hmm_suffix + " > tmpOrAnPipeHMMs/" + file_name + hmm_suffix + ".lenAdded.hmm", shell = True)
    if not args.thammerin_results and (args.search_mode == "fl" or not "fullLenForHMM" in file_name):
        if not "tmpOrAnPipeFrames.fa" in os.listdir('./') and not args.genome_orfs:
            subprocess.call(thammerin + " --frames_out tmpOrAnPipeFrames.fa -e " +
                             str(args.evalue) + " -p tmpOrAnPipeHMMs/" +
                             file_name + hmm_suffix + ".lenAdded.hmm -n " + target_genome +
                             ' > ' + file_name[:-3] + '.tmpOrAnPipe.thammerin.tab',
                             shell = True)
        else:
            subprocess.call(thammerin + " --frames_in " + thammerin_orfs + " -e " + str(args.evalue) + " -p tmpOrAnPipeHMMs/" + file_name + hmm_suffix + '.lenAdded.hmm > ' + file_name[:-3] + '.tmpOrAnPipe.thammerin.tab', shell = True)
    #debug
        print "Running " + thammerin + " --frames_in " + thammerin_orfs + " -e " + str(args.evalue) + \
            " -p tmpOrAnPipeHMMs/" + file_name + hmm_suffix + '.lenAdded.hmm'
    #/debug

#Writes protein sequences, builds a distance tree using mafft, and parses this tree to give hmm groups
if (dthresh < 1 and args.protein_seqs) or args.annotations:
    clusters = build_clusters()
    output_fastas(args.search_mode)
    query_pep_files = os.listdir('tmpOrAnPipeHMMs')
    query_pep_dir = 'tmpOrAnPipeHMMs/'
elif dthresh < 1:
    sys.stderr.write("Warning: --cutoff was less than 1, indicating that you wanted HAP to split queries into clusters, \
    however query sequences were not provided in a splitable format. Continuing assuming --cutoff value was an error.\n")
#No clustering needed, proceeds
elif args.protein_seqs:
    query_pep_files = [protein_seqs.split('/')[-1]]
    query_pep_dir = "/".join(protein_seqs.split('/')[:-1]) + '/'
elif args.hmm:
    query_pep_files = [protein_hmm.split('/')[-1]]
    query_pep_dir = "/".join(protein_hmm.split('/')[:-1]) + '/'
    hmm_dir = ""
elif args.fasta_dir:
    query_pep_files = os.listdir(fasta_dir)
    query_pep_dir = fasta_dir
elif args.alignment_dir:
    query_pep_files = os.listdir(alignment_dir)
    query_pep_dir = alignment_dir
elif args.hmm_dir:
    query_pep_files = os.listdir(hmm_dir)
    query_pep_dir = hmm_dir



for file_name in query_pep_files:
    if not args.hmm_dir and not args.hmm:
        hmm_dir = "tmpOrAnPipeHMMs/"
    if not "tmpOrAnPipeFrames.fa" in os.listdir('./') and not args.genome_orfs:
        thammerin_func(file_name,hmm_dir)
    else:
        threading.Thread(target = thammerin_func,args = [file_name,hmm_dir]).start()
        while threading.active_count() >= 2 + args.threads / 3:
            time.sleep(0.1)

if not args.thammerin_results:
    while threading.active_count() > 1:
        time.sleep(1)
    subprocess.call('cat *.tmpOrAnPipe.thammerin.tab > tmpOrAnPipe.thammerin.tab', shell = True)
    blasttab = 'tmpOrAnPipe.thammerin.tab'
else:
    blasttab = call_directory + '/' + args.thammerin_results

#Collects thammerin results and parses to generate candidate loci
coord_dict = {}

blast_tab = open(blasttab).read().split('\n')
if blast_tab[-1] == "":
    blast_tab.pop()

if args.search_mode == "exons":
    coord_dict = CandidateLociBuilder.build_coord_dict_based_on_exon_number(blast_tab)
elif args.search_mode == "fl":
    coord_dict = CandidateLociBuilder.build_coord_dict_based_on_seq_length(blast_tab)

candidate_loci_fasta_list = CandidateLociBuilder.find_candidate_loci(coord_dict, args.buffer, target_genome)

candidate_loci = open('candidate_loci.fasta','w')
candidate_loci.write('\n'.join(candidate_loci_fasta_list))
candidate_loci.close()


if args.annotator == "ABCENTH":
    subprocess.call('python ' + program_dir +
                    'HitTabFilter.py --tab tmpOrAnPipe.thammerin.tab > thammerin.unique.tab', shell = True)
    subprocess.call('python ' + program_dir +
                    'ParseHAPpyTable.py --tab thammerin.unique.tab > thammerin.unique_forABCENTH.tab', shell = True)
    subprocess.call('python ' + program_dir +
                    'ParseHAPpyTable.py --hmm_dir tmpOrAnPipeHMMs/ > clusterInfo.tab', shell = True)
    subprocess.call('python ' + program_dir +
                    'ABCENTH.py --table thammerin.unique_forABCENTH.tab --query_exon_info_table clusterInfo.tab --genome ' +
                    target_genome + ' --log ABCENTH.log > ../ABCENTH.gtf',shell = True)
else:
    candidate_loci_dict = {}
    if args.search_mode == 'exons':
        file_key = "fullLenForHMM.fa.hmm.lenAdded.hmm"
    else:
        file_key = "lenAdded.hmm"
    for locus in candidate_loci_fasta_list:
        candidate_loci_dict[locus.split('\n')[0][1:]] = locus.split('\n')[1]
    AnnotatorRunner.parse_loci('candidate_loci.fasta',candidate_loci_dict,"tmpOrAnPipeHMMs/","Clusters",
                              thammerin,genewisedb, file_key = file_key, threads = args.threads,
                              annotator = args.annotator,augustus_profile_dir = augustus_profile_dir,
                              augustus_cfg = "--extrinsicCfgFile=" + program_dir + 'augustus_extrinsic.MP.cfg',
                              augustus_species = args.augustus_species)
    if 'augustus' in args.annotator:
        subprocess.call('cat Clusters/*/augustus.gtf > HAPpy_augustus.gtf',shell = True)
        subprocess.call('cat Clusters/*/augustus.pep > HAPpy_augustus.pep',shell = True)
    if 'genewise' in args.annotator:
        subprocess.call('cat Clusters/*/genewise.gtf > HAPpy_genewise.gtf',shell = True)
        subprocess.call('cat Clusters/*/genewise.pep > HAPpy_genewise.pep',shell = True)
