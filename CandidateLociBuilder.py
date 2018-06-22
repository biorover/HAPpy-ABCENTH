#!/usr/bin/env python

#Library for parsing homology hit tables to find candidate loci

import copy

def build_coord_dict_based_on_exon_number(hit_table):
    """finds any hits in a table of homology hits that seem to be first, last, or almost first/last exons, based on
    exon numbers given in name"""
    coord_dict = {}
    for hsp in hit_table:
        fields = hsp.split('\t')
        tid = fields[1]
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
    return coord_dict


def build_coord_dict_based_on_seq_length(hit_table, start_thresh = 10, almost_start_thresh = 40,
                                         end_thresh = 5, almost_end_thresh = 30):
    """finds any hits in a table of homology hits that seem to be first, last, or almost first/last exons, based on
    hit coordinates and the sequence length given in the name ("_len_NN")"""
    coord_dict = {}
    for hsp in hit_table:
        fields = hsp.split('\t')
        tid = fields[1]
        qlen = int(fields[0].split('_len_')[1])
        dict_add = []
        if int(fields[6]) < start_thresh:
            if int(fields[9]) > int(fields[8]):
                dict_add.append((int(fields[8]), "start"))
            else:
                dict_add.append((int(fields[8]), "end"))
        elif int(fields[6]) < almost_start_thresh:
            if int(fields[9]) > int(fields[8]):
                dict_add.append((int(fields[8]), "almost_start"))
            else:
                dict_add.append((int(fields[8]), "almost_end"))
        if qlen - int(fields[7]) < end_thresh:
            if int(fields[9]) > int(fields[8]):
                dict_add.append((int(fields[9]), "end"))
            else:
                dict_add.append((int(fields[9]), "start"))
        elif qlen - int(fields[7]) < almost_end_thresh:
            if int(fields[9]) > int(fields[8]):
                dict_add.append((int(fields[9]), "almost_end"))
            else:
                dict_add.append((int(fields[9]), "almost_start"))
        if dict_add != []:
            if tid in coord_dict:
                coord_dict[tid].extend(copy.deepcopy(dict_add))
            else:
                coord_dict[tid] = copy.deepcopy(dict_add)
    return coord_dict

def find_candidate_loci(coord_dict, buffer_len, genome_file):
    """break into candidate loci (between start HSP and end HSP, +- buffer bp)"""
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
                elif coord_dict[tid][coord_index][1] == "almost_start" and (coord_dict[tid][coord_index - 1][1] in ['end','almost_end'] or coord_index == 0):
                    counter = coord_index + 1
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
                if len(candidate_loci[tid]) > 0:
                    largest_end = candidate_loci[tid][-1][1]

    genome_dict = {}
    genome_list = open(genome_file).read().split('>')[1:]
    for seq in genome_list:
        lines = seq.split('\n')
        genome_dict[lines[0]] = "".join(lines[1:])

    candidate_loci_seqs = {}
    candidate_loci_fasta_list = []
    for tid in candidate_loci:
        scaf_seqlen = len(genome_dict[tid])
        for locus in candidate_loci[tid]:
            if locus[0] > buffer_len and buffer_len + locus[1] < scaf_seqlen:
                locus_seq = genome_dict[tid][locus[0] - buffer_len : locus[1] + buffer_len]
                locus_start = locus[0] - buffer_len
                locus_stop = locus[1] + buffer_len
            elif locus[0] <= buffer_len and buffer_len + locus[1] < scaf_seqlen:
                locus_seq = genome_dict[tid][: locus[1] + buffer_len]
                locus_start = 0
                locus_stop = locus[1] + buffer_len
            elif buffer_len + locus[1] >= scaf_seqlen and locus[0] > buffer_len:
                locus_seq = genome_dict[tid][locus[0] - buffer_len :]
                locus_start = locus[0] - buffer_len
                locus_stop = scaf_seqlen
            else:
                locus_seq = genome_dict[tid][:]
                locus_start = 0
                locus_stop = scaf_seqlen
            candidate_loci_seqs[tid + "_" + str(locus[0]) + '-' + str(locus[1])] = locus_seq
            candidate_loci_fasta_list.append('>' + tid + "_" + str(locus_start) + '-' + str(locus_stop) + '\n' + locus_seq + '\n')
    return candidate_loci_fasta_list