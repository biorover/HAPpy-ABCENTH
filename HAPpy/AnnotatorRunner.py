#!/usr/bin/env python

#Seperates candidate loci by cluster and annotates them with genewise

import subprocess
import os
import tempfile
from HAPpy import genome_fork as genome
from HAPpy import toolbox_for_HAP
import threading
import time

def thammerin_for_parse_loci(hmm_file,hmm_dir,thammerin_path,candidate_loci_orfs_name):
    subprocess.call(thammerin_path + " -p " + hmm_dir + "/" + hmm_file +
                              " --frames_in " + candidate_loci_orfs_name + ' > thammerin_canloc_' +
                              hmm_file + '.tab', shell = True)


def parse_loci(candidate_loci_file,candidate_loci_seqs,hmm_dir,output_dir,thammerin_path = "thammerin.py",
               genewisedb_path = "genewisedb",hmmconvert_path = "hmmconvert", augustus_path = 'augustus',
               augustus_profile_dir = None, file_key = ".hmm", annotator = 'genewise', threads = 1, augustus_cfg = " ",
               augustus_species = 'fly'):
    candidate_loci_orfs = tempfile.NamedTemporaryFile()
    orfs_file_written = False
    for hmm_file in os.listdir(hmm_dir):
        if file_key in hmm_file:
            subprocess.call("mkdir -p " + output_dir + '/' + hmm_file[:-4], shell = True)
            if not orfs_file_written:
                subprocess.call(thammerin_path + " -p " + hmm_dir + "/" + hmm_file +
                                                  " -n " + candidate_loci_file + " --frames_out "
                                                  + candidate_loci_orfs.name + ' > thammerin_canloc_' + hmm_file + '.tab',
                                                  shell = True)
                orfs_file_written = True
            else:
                threading.Thread(target = thammerin_for_parse_loci, args =
                                 [hmm_file,hmm_dir,thammerin_path,candidate_loci_orfs.name]).start()
                while threading.active_count() > 1 + threads / 3:
                    time.sleep(1)
    while threading.active_count() > 1:
        time.sleep(1)
    subprocess.call('cat thammerin_canloc_* > thammerin_allcanloc.tab', shell = True)
    thammerin_lines = open('thammerin_allcanloc.tab').read().split('\n')
    candidate_loci_orfs.close()
    while "" in thammerin_lines:
        thammerin_lines.remove("")
    thammerin_scores = {}
    thammerin_hit_dict = {}
    for hit in thammerin_lines:
        fields = hit.split('\t')
        class_hit = fields[0].split('_hmmFile_')[1].split('_len_')[0][:-4]
        loc_name = fields[1]
        score = float(fields[11])
        if loc_name in thammerin_scores:
            if score > thammerin_scores[loc_name][0]:
                thammerin_scores[loc_name] = [score,class_hit]
            if class_hit in thammerin_hit_dict[loc_name]:
                thammerin_hit_dict[loc_name][class_hit].append(hit)
            else:
                thammerin_hit_dict[loc_name][class_hit] = [hit]
        else:
            thammerin_hit_dict[loc_name] = {class_hit:[hit]}
            thammerin_scores[loc_name] = [score,class_hit]
    for locus in thammerin_scores:
        cluster = thammerin_scores[locus][1]
        class_fasta = open(output_dir + "/" + cluster + '/candidate_loci.fasta','a')
        class_fasta.write('>' + locus + '\n' + candidate_loci_seqs[locus] + '\n')
        class_fasta.close()
        #locus_fasta = open(output_dir + "/" + cluster + '/' + locus + ".fa",'w')
        #locus_fasta.write('>' + locus + '\n' + candidate_loci_seqs[locus] + '\n')
        #locus_fasta.close()
        #starts prepping hints
    for cluster_folder in os.listdir(output_dir):
        threading.Thread(target= annotate, args = [annotator,hmm_dir,cluster_folder,output_dir,
                                                   thammerin_hit_dict,genewisedb_path,augustus_path,
                                                   hmmconvert_path, augustus_profile_dir, augustus_cfg,
                                                   augustus_species]).start()
        while threading.active_count() >= 2 + threads:
            time.sleep(1)



def annotate(annotator,hmm_dir,cluster_folder,output_dir,thammerin_hit_dict,genewisedb_path,
             augustus_path,hmmconvert_path,augustus_profile_dir, augustus_cfg, augustus_species):
    hints_file = open(output_dir + "/" + cluster_folder+ '/augustus_hints.gff','a')
    loci_list = subprocess.check_output('grep ">" ' + output_dir + '/' + cluster_folder + '/candidate_loci.fasta',
                                        shell = True).split('\n')
    while "" in loci_list:
        loci_list.remove("")
    for raw_locus in loci_list:
        locus = raw_locus[1:]
        hit_annotation = genome.read_blast_csv("\n".join(thammerin_hit_dict[locus][cluster_folder]).replace('\t',','))
        hints = genome.write_gff(hit_annotation,'augustus hint exonpart P').split('\n')
        for hint_line_index in range(len(hints)):
            hint_line_fields = hints[hint_line_index].split('\t')
            if int(hint_line_fields[4]) - int(hint_line_fields[3]) > 20:
                hint_line_fields[3] = str(int(hint_line_fields[3]) + 10)
                hint_line_fields[4] = str(int(hint_line_fields[4]) - 10)
                hints[hint_line_index] = "\t".join(hint_line_fields)
        hints_file.write('\n'.join(hints) + '\n')
    subprocess.call(hmmconvert_path + ' -2 ' + hmm_dir + '/' + cluster_folder + ".hmm > " + output_dir + '/' +
                   cluster_folder + '/' + cluster_folder + '.hmm', shell = True)
    if 'genewise' in annotator:
        subprocess.call(genewisedb_path + ' -cut 10 -pseudo -pep -gff -hmmer ' + output_dir + '/' +
                        cluster_folder + '/' + cluster_folder + '.hmm ' +
                        output_dir + '/' + cluster_folder + '/candidate_loci.fasta > ' + output_dir +
                        '/' + cluster_folder + '/candidate_loci.genewise 2>> genewise_log.txt' ,shell = True)
        genewise_lines,genewise_fasta = toolbox_for_HAP.parse_genewise(output_dir + "/" + cluster_folder +
                                                                       '/candidate_loci.genewise',
                                                                       name_prefix = cluster_folder.replace('.lenAdded','').replace('.hmm',''))
        ggtf_file = open(output_dir + "/" + cluster_folder + '/genewise.gtf','w')
        ggtf_file.write(''.join(genewise_lines))
        ggtf_file.close()
        gfasta_file = open(output_dir + "/" + cluster_folder + '/genewise.pep','w')
        gfasta_file.write(''.join(genewise_fasta))
        gfasta_file.close()
        if 'augustus' in annotator:
            genewiselines = open(output_dir + '/' + cluster_folder + '/candidate_loci.genewise').read().split('\n')
            hints_file.write('\n'.join(toolbox_for_HAP.genewise2aughints(genewiselines)) + '\n')
    hints_file.close()
    if 'augustus' in annotator:
        if augustus_profile_dir:
            augprofstring = " --proteinprofile=" + augustus_profile_dir + '/' + cluster_folder.replace('.hmm.lenAdded','') + '.prfl '
        else:
            augprofstring = " "
        print("running " + augustus_path + augprofstring + augustus_cfg + " --species=" + augustus_species + \
                        " --hintsfile=" + output_dir + "/" + cluster_folder + \
                        '/augustus_hints.gff '+ output_dir + '/' + cluster_folder + '/candidate_loci.fasta')
        subprocess.call(augustus_path + augprofstring + augustus_cfg + " --species=" + augustus_species +
                        " --hintsfile=" + output_dir + "/" + cluster_folder +
                        '/augustus_hints.gff '+ output_dir + '/' + cluster_folder + '/candidate_loci.fasta >> ' +
                        output_dir + "/" + cluster_folder + '/augustus.out', shell = True)
        augustus_lines,augustus_fasta = toolbox_for_HAP.parse_augustus(output_dir + "/" + cluster_folder +
                                                                       '/augustus.out', name_prefix = cluster_folder.replace('lenAdded','').replace('.hmm',''))
        gtf_file = open(output_dir + "/" + cluster_folder + '/augustus.gtf','w')
        gtf_file.write(''.join(augustus_lines))
        gtf_file.close()
        fasta_file = open(output_dir + "/" + cluster_folder + '/augustus.pep','w')
        fasta_file.write(''.join(augustus_fasta))
        fasta_file.close()
