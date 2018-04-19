#!/usr/bin/env python

import sys
import os
import argparse
import genome_fork as genome
import subprocess
import ete2
import time
import copy

parser = argparse.ArgumentParser(description="Pipeline for annotating ORs in a genome using homology and exon \
                                 phase information from a closely related species.\n\nProgram dependencies:\n\
                                 python\nmafft\nhmmer suite\nthammerin.py\nGeMoMa.\n\nAdditionally, this program requires \
                                 the genome library from MAGOT (multitool for analyzing genomes or \
                                 transcriptomes, github.com/biorover/MAGOT) and the ete2 python library. These \
                                 libraries must be in a directory in the PYTHONPATH \
                                 environmental variable. If you don't know what that means, just install ete2 \
                                 according to the directions on the programs website, and download MAGOT \
                                 to some directory and in any given terminal session in which you'll use this\
                                 tool, type \"export PYTHONPATH=$PYTHONPATH:/path/to/MAGOTdirectory\", or add \
                                 that command to your .profile (being sure to source your profile after you \
                                 add it). One day I promise to set up MAGOT with an install script. Maybe by the \
                                 time you're using this I'll have done that! Wow, you have a high tolerance for \
                                 reading long help texts. Good job!")

parser.add_argument('-a','--annotations',dest='annotations', 
                    help = 'Annotations (in gtf format)')

parser.add_argument('-r','--ref_genome',dest="ref_genome",
                    help = 'genome of closely related species')

parser.add_argument('-g','--genome',dest='target_genome', 
                    help = 'genome to be annotated')

parser.add_argument('-f','--program_filepaths',dest = 'program_filepaths', default = None, 
                    help = 'optional list of file paths for programs not in PATH variable, formated as "programName=/path/to/programName"')

parser.add_argument('-m','--min_orf_size',dest='min_orf_size', default = 10,
                    help = 'minimum size for orfs to be search by hmmsearch')

parser.add_argument('-c','--cutoff',dest = 'distance_cutoff', default = 0.5, type = float,
                    help = 'distance cutoff for seperating proteins into clusters for hmm building. Will have\
                           to be determined empirically.')

parser.add_argument('-t','--threads', dest = 'threads', default = 1, type = int,
                    help = "number of threads to be used with processes that support multithreading (currently just mafft)")

parser.add_argument('--tree',dest = 'gene_tree', default = None,
                    help = 'pre-computed UPGMA gene tree for query genes to save time')

parser.add_argument('--thammerin_results', dest = 'thammerin_results', default = None,
                    help = 'pre-computed thammerin results')

parser.add_argument('--buffer', dest = 'buffer', default = 500,
                    help = 'buffer on either side of loci identified to feed into gene predictor')

parser.add_argument('--evalue', default = 0.05, type = float, help = 'Evalue cutoff for thammerin')

parser.add_argument('--genome_orfs', default = None, help = 'ORFs file from previous thammerin run (saves about five minutes for insect-sized genomes)')

args = parser.parse_args()


#sets up some basic variables
dthresh = args.distance_cutoff
if "genome_orfs" in args:
    thammerin_orfs = os.path.abspath(args.genome_orfs)
else:
    thammerin_orfs = 'tmpOrAnPipeFrames.fa'

#sets up default program paths, overwritten by any program paths passes with -f or --program_filepaths
hmmsearch = "hmmsearch"
hmmbuild = "hmmbuild"
mafft = "mafft"
thammerin = "thammerin.py"

if args.program_filepaths:
    for line in args.program_filepaths:
        exec(line.replace('\n','').replace('\r',''))


#Reads in annotations, sets annotation genome path
ref_genomes = []
genome_paths = args.ref_genome.split(',')
geneset_paths = args.annotations.split(',')
for genome_index in range(len(genome_paths)):
    ref_genomes.append(genome.Genome(genome_paths[genome_index]))
    ref_genomes[genome_index].read_gff(geneset_paths[genome_index])

    
target_genome = os.path.abspath(args.target_genome)
call_directory = os.getcwd()

subprocess.call('mkdir -p tmpOrAnPipeDir',shell = True)
os.chdir('tmpOrAnPipeDir')
#Writes protein sequences, builds a distance tree using mafft, and parses this tree to give hmm groups
flprots = open('tmpOrAnPipe.flpeps.fa','w')
for ref_genome in ref_genomes:
    flprots.write(ref_genome.annotations.get_fasta('gene',seq_type = 'protein')+'\n')

flprots.close()

if not args.gene_tree:
    subprocess.call(mafft + ' --thread ' + str(args.threads) + ' --globalpair --treeout tmpOrAnPipe.flpeps.fa > tmpeOrAnPipe.tmp', shell = True)
    prottree = ete2.Tree(open('tmpOrAnPipe.flpeps.fa.tree').read() + ';')
else:
    prottree = ete2.Tree(open(call_directory + '/' + args.gene_tree).read())
    
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

#Builds HMMs for exons in clusters and runs thammerin with all hmms vs the genome to be annotated
cluster_dict = {}
lengths_dict = {}
subprocess.call('mkdir tmpOrAnPipeHMMs', shell=True)
for cluster_index in range(len(clusters)):
    cluster = clusters[cluster_index]
    cluster_dict[cluster_index] = []
    exon_number = None
    exon_phases = None
    exon_lengths = None
    exon_file_root = "tmpOrAnPipeHMMs/cluster" + str(cluster_index)
    for seqname in cluster:   
        transcript_id = "_".join(seqname.split("_")[1:])
        for ref_genome in ref_genomes:
            try:
                cluster_dict[cluster_index].append(ref_genome.annotations.transcript[transcript_id])
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
                lengths_dict[cluster_index] = exon_lengths[:]
            except KeyError:
                pass

running_commands = [] #for multithreading

exon_pep_files = os.listdir("tmpOrAnPipeHMMs/")

for file_name in exon_pep_files:
    if open("tmpOrAnPipeHMMs/" + file_name).read().count('>') > 1:
        print "Aligning " + file_name
        subprocess.call("mafft --globalpair --maxiterate 1000 tmpOrAnPipeHMMs/" + file_name + " > tmpOrAnPipeHMMs/"+ file_name[:-3] + ".mafftGinsi.fa 2> tmp.tmp" , shell = True)
        hmmer_input = file_name[:-3] + ".mafftGinsi.fa"
    else:
        print "Passing along " + file_name
        hmmer_input = file_name
    try:
        query_len = subprocess.check_output(hmmbuild + " --amino tmpOrAnPipeHMMs/" + file_name[:-3] + ".hmm tmpOrAnPipeHMMs/" +  hmmer_input, shell = True).split('\n')[12].split()[4]
    except subprocess.CalledProcessError:
        print "Hm, problem with hmmbuild for file " + file_name
    gene_cluster = file_name.split('exon')[0]
    subprocess.call("sed -i -e 's/mafftGinsi/___" + gene_cluster + "_len_" + query_len + "/g' tmpOrAnPipeHMMs/" + file_name[:-3] + ".hmm", shell = True)
    
    if not args.thammerin_results:
        if not "tmpOrAnPipeFrames.fa" in os.listdir('./') and not "genome_orfs" in args:
            running_commands.append(subprocess.Popen(thammerin + " --frames_out tmpOrAnPipeFrames.fa -e " + str(args.evalue) + " -p tmpOrAnPipeHMMs/" + file_name[:-3] + ".hmm -n " + target_genome + ' > ' + file_name[:-3] + '.tmpOrAnPipe.thammerin.tab', shell = True))
            for command in running_commands:
                command.wait()
        else:
            running_commands.append(subprocess.Popen(thammerin + " --frames_in " + thammerin_orfs + " -e " + str(args.evalue) + " -p tmpOrAnPipeHMMs/" + file_name[:-3] + '.hmm > ' + file_name[:-3] + '.tmpOrAnPipe.thammerin.tab', shell = True))
    #debug
        print "Running " + thammerin + " -e " + str(args.evalue) + " -p tmpOrAnPipeHMMs/" + file_name[:-3] + ".hmm -n " + target_genome
    #/debug
    ongoing_commands_count = 0
    for cmd in running_commands:
        if cmd.poll() == None:
            ongoing_commands_count = ongoing_commands_count + 1
    while ongoing_commands_count >= args.threads / 3:
        time.sleep(5)
        ongoing_commands_count = 0
        for cmd in running_commands:
            if cmd.poll() == None:
                ongoing_commands_count = ongoing_commands_count + 1

if not args.thammerin_results:
    for command in running_commands:
        command.wait()

    subprocess.call('cat *.tmpOrAnPipe.thammerin.tab > tmpOrAnPipe.thammerin.tab', shell = True)
    blasttab = 'tmpOrAnPipe.thammerin.tab'
else:
    blasttab = call_directory + '/' + args.thammerin_results

#Collects thammerin results and parses to generate candidate loci
coord_dict = {}

blast_tab = open(blasttab).read().split('\n')
if blast_tab[-1] == "":
    blast_tab.pop()

for hsp in blast_tab:
    fields = hsp.split('\t')
    tid = fields[1]
    #qlen = int(fields[0].split('_len_')[1])
    dist_from_start = int(fields[0].split('exon')[1].split('of')[0])
    dist_from_end = int(fields[0].split('of')[1].split('phase')[0]) - dist_from_start
    dict_add = []
    if dist_from_start == 0:
        if int(fields[9]) > int(fields[8]):
            dict_add.append((int(fields[8]), "start"))
        else:
            dict_add.append((int(fields[8]), "end"))
    elif dist_from_start == 1:
        if int(fields[9]) > int(fields[8]):
            dict_add.append((int(fields[8]), "almost_start"))
        else:
            dict_add.append((int(fields[8]), "almost_end"))
    if dist_from_end == 0:
        if int(fields[9]) > int(fields[8]):
            dict_add.append((int(fields[9]), "end"))
        else:
            dict_add.append((int(fields[9]), "start"))
    elif dist_from_end == 1:
        if int(fields[9]) > int(fields[8]):
            dict_add.append((int(fields[9]), "almost_end"))
        else:
            dict_add.append((int(fields[9]), "almost_start"))
    if dict_add != []:
        if tid in coord_dict:
            coord_dict[tid].extend(copy.deepcopy(dict_add))
        else:
            coord_dict[tid] = copy.deepcopy(dict_add)


#break into candidate loci (between start HSP and end HSP, +- buffer bp)
candidate_hits = open("Candidate_start-stop_hits.txt",'w')
for tid in coord_dict:
    coord_dict[tid].sort()
    candidate_hits.write(tid + ': ' + str(coord_dict[tid]) + '\n')
candidate_hits.close()

candidate_loci = {}
for tid in coord_dict:
    candidate_loci[tid] = []
    coord_dict[tid].sort()
    largest_end = 0
    for coord_index in range(len(coord_dict[tid]) - 1):
        if coord_dict[tid][coord_index + 1][0] - coord_dict[tid][coord_index][0] < 10000 and coord_dict[tid][coord_index][0] > largest_end:
            if coord_dict[tid][coord_index][1] == "start" and coord_dict[tid][coord_index + 1][1] == 'end':
                candidate_loci[tid].append((coord_dict[tid][coord_index][0],coord_dict[tid][coord_index + 1][0]))
            elif coord_dict[tid][coord_index][1] == "start" and coord_dict[tid][coord_index + 1][1] == 'almost_start' and coord_index + 2 < len(coord_dict[tid]):
                counter = coord_index + 2
                end_found = False
                almost_end_found = False
                while coord_dict[tid][counter][0] - coord_dict[tid][coord_index][0] < 10000 and not coord_dict[tid][counter][1] == "start":
                    if coord_dict[tid][counter][1] == 'end':
                        candidate_loci[tid].append((coord_dict[tid][coord_index][0],coord_dict[tid][counter][0]))
                        end_found = True
                        break
                    elif counter + 1 == len(coord_dict[tid]):
                        break
                    elif coord_dict[tid][counter][1] == 'almost_end':
                        highest_almost_end_index = counter
                        almost_end_found = True
                    counter = counter + 1
                if almost_end_found and not end_found:
                    candidate_loci[tid].append((coord_dict[tid][coord_index][0],coord_dict[tid][highest_almost_end_index][0]))
            elif coord_dict[tid][coord_index][1] == "start" and coord_dict[tid][coord_index + 1][1] == 'almost_end' and coord_index + 2 < len(coord_dict[tid]):
                counter = coord_index + 2
                highest_almost_end_index = coord_index + 1
                end_found = False
                while coord_dict[tid][counter][0] - coord_dict[tid][coord_index][0] < 10000 and not coord_dict[tid][counter][1] in ["start","almost_start"]:
                    if coord_dict[tid][counter][1] == 'end':
                        candidate_loci[tid].append((coord_dict[tid][coord_index][0],coord_dict[tid][counter][0]))
                        end_found = True
                        break
                    elif counter + 1 == len(coord_dict[tid]):
                        break
                    elif coord_dict[tid][counter][1] == 'almost_end':
                        highest_almost_end_index = counter
                    counter = counter + 1
                if not end_found:
                    candidate_loci[tid].append((coord_dict[tid][coord_index][0],coord_dict[tid][highest_almost_end_index][0]))
            elif coord_dict[tid][coord_index][1] == "almost_start" and coord_dict[tid][coord_index - 1][1] in ['end','almost_end']:
                counter = coord_index + 1
                end_found = False
                almost_end_found = False
                while coord_dict[tid][counter][0] - coord_dict[tid][coord_index][0] < 10000 and not coord_dict[tid][counter][1] in ["start","almost_start"]:
                    if coord_dict[tid][counter][1] == 'end':
                        candidate_loci[tid].append((coord_dict[tid][coord_index][0],coord_dict[tid][counter][0]))
                        end_found = True
                        break
                    elif counter + 1 == len(coord_dict[tid]):
                        break
                    elif coord_dict[tid][counter][1] == 'almost_end':
                        highest_almost_end_index = counter
                        almost_end_found = True
                    counter = counter + 1
                if almost_end_found and not end_found:
                    candidate_loci[tid].append((coord_dict[tid][coord_index][0],coord_dict[tid][highest_almost_end_index][0]))
            if len(candidate_loci[tid]) > 0:
                largest_end = candidate_loci[tid][-1][1]

genome_dict = {}
genome_list = open(target_genome).read().split('>')[1:]
for seq in genome_list:
    lines = seq.split('\n')
    genome_dict[lines[0]] = "".join(lines[1:])


candidate_loci_seqs = {}
candidate_loci_fasta = open('candidate_loci.fasta','w')
for tid in candidate_loci:
    scaf_seqlen = len(genome_dict[tid])
    for locus in candidate_loci[tid]:
        if locus[0] > args.buffer and args.buffer + locus[1] < scaf_seqlen:
            locus_seq = genome_dict[tid][locus[0] - args.buffer : locus[1] + args.buffer]
        elif locus[0] <= args.buffer and args.buffer + locus[1] < scaf_seqlen:
            locus_seq = genome_dict[tid][: locus[1] + args.buffer]
        elif args.buffer + locus[1] >= scaf_seqlen and locus[0] > args.buffer:
            locus_seq = genome_dict[tid][locus[0] - args.buffer :]
        else:
            locus_seq = genome_dict[tid][:]
        candidate_loci_seqs[tid + "_" + str(locus[0]) + '-' + str(locus[1])] = locus_seq
        candidate_loci_fasta.write('>' + tid + "_" + str(locus[0]) + '-' + str(locus[1]) + '\n' + locus_seq + '\n')
candidate_loci_fasta.close()




#These parameters will be needed for GeMoMa tblastn: tblastn -outfmt "6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles"