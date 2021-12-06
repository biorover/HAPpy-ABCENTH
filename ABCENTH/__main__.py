#!/usr/bin/env python

#   Features to be added later:
#proper splice site modeler (perhaps optimizing splice site score vs. difference from expected length), maybe MAYBE taking
#into account non-conserved exon phase.



import argparse
from HAPpy import genome_fork as genome
from ABCENTH.FindOrfs import orf_finder,genewisesearch
import sys, os

def isyes(x):
    if x in [True,'True','true','TRUE','TrUe','T','t','oui','si','y','Y','Yes','YES','yes','oh my YESSSSS!']:
        return True
    else:
        return False


#The heart of gold of this program- the function to define an exon based on the coords and target phase
def exon_finder(tstart,tend,strand,qstart,qend,qlen,qstartphase,qendphase,seqdict,seqname,
                max_offset = 30, is_start = False, is_stop = False, nevermind_atg = False,
                cluster = None,exon_number = None,log_file = open(os.devnull, 'w'),full_pseudoexon_search = True,
                exon_info_dict = False):
    """"finds exons with ORFs based on the requested parameters. Note that it is expected that tstart < tend, \
    even for minus strand features, so these might need to be reversed if coming from say tblastn. \
    Also coords are expected as 1-based (as output from blast), and are \
    converted internally to 0 based. It is also expected that the hit itself doesn't contain stop codons."""
    start = None
    end = None
    pseudo = False
    max_coord = len(seqdict[seqname]) - 1
    if strand == "+":
        phasestart_offset = (3 - qstartphase) % 3
        phasestop_offset = qendphase
        start_match_offset, stop_match_offset = 3 * (qstart - 1),3 * (qlen - qend)
    elif strand == '-':
        phasestart_offset = qendphase
        phasestop_offset = (3 - qstartphase) % 3
        start_match_offset, stop_match_offset = 3 * (qlen - qend), 3 * qstart
    ideal_start = tstart - phasestart_offset - start_match_offset
    ideal_end = tend + stop_match_offset + phasestop_offset
    pseudo_start = tstart - phasestart_offset
    pseudo_end = tend + phasestop_offset
    gc_start, gc_end = None, None
    for offset in range(0,max_offset + 3,3):
        if start:
            break
        for direction in [1,-1]:
            test_start = ideal_start - offset * direction
            test_seq = genome.Sequence(seqdict[seqname][test_start-1 + phasestart_offset:tend])
            if strand == "-":
                test_seq = test_seq.reverse_compliment()
            if not test_seq.translate():
                continue
            elif is_stop and strand == "-":
                if ideal_start - 1 < test_start < pseudo_start:
                    pseudo_start = test_start
                lastcodon = seqdict[seqname][test_start - 1:test_start + 2]
                if lastcodon.upper() in ['TTA','TCA','CTA'] and test_seq.translate().count('*') == 1:
                    start = test_start
                    break
            elif not "*" in test_seq.translate():
                if ideal_start - 1 < test_start < pseudo_start:
                    pseudo_start = test_start
                if is_start and strand == '+':
                    if nevermind_atg:
                        start = test_start
                        break
                    else:
                        firstcodon = seqdict[seqname][test_start - 1:test_start + 2]
                        if firstcodon.upper() == "ATG":
                            start = test_start
                            break
                else:
                    splicesite = seqdict[seqname][test_start-3:test_start-1]
                    if (strand == '+' and splicesite.upper() == "AG") or (strand == '-' and splicesite.upper() == "AC"):
                        start = test_start
                        break
                    elif strand == '-' and splicesite.upper() == "GC" and not gc_start:
                        gc_start = test_start
    if not start:
        if gc_start:
            start = gc_start
        else:
            pseudo = "P"
            start = pseudo_start
    for offset in range(0,max_offset + 3,3):
        if end:
            break
        for direction in [-1,1]:
            test_end = ideal_end - offset * direction
            if test_end - start < 3:
                break
            test_seq = genome.Sequence(seqdict[seqname][start - 1 + phasestart_offset:test_end - phasestop_offset])
            if strand == "-":
                test_seq = test_seq.reverse_compliment()
            if not test_seq.translate():
                continue
            elif is_stop and strand == "+":
                if ideal_end + 1 > test_end > pseudo_end:
                    pseudo_end = test_end
                lastcodon = seqdict[seqname][test_end - 3:test_end]
                if lastcodon.upper() in ['TAA','TGA','TAG'] and test_seq.translate().count('*') == 1:
                    end = test_end
                    break
            elif not "*" in test_seq.translate() or (is_stop and not "*" in test_seq.translate()[:-1]):
                if ideal_end + 1 > test_end > pseudo_end:
                    pseudo_end = test_end
                if is_start and strand == '-':
                    if nevermind_atg:
                        end = test_end
                        break
                    else:
                        firstcodon = seqdict[seqname][test_end - 3:test_end]
                        if firstcodon.upper() == "CAT":
                            end = test_end
                            break
                else:
                    splicesite = seqdict[seqname][test_end:test_end + 2]
                    if (strand == '+' and splicesite.upper() == "GT") or (strand == '-' and splicesite.upper() == "CT"):
                        end = test_end
                        break
                    elif strand == "+" and splicesite.upper() == "GC" and not gc_end:
                        gc_end = test_end
    if not end:
        if gc_end:
            end = gc_end
        else:
            pseudo = "P"
            end = pseudo_end
    start = max([1,start])
    end = min([end,max_coord])
    if pseudo and full_pseudoexon_search and cluster != None and exon_number != None and exon_info_dict:
        gwexons = genewisesearch(seqdict[seqname],qstartphase,qendphase,strand,
                    exon_info_dict[str(cluster) + ':' + str(exon_number)][6], 
                    search_coords = [ideal_start - 3 - max_offset,ideal_end + 3 + max_offset],
                    seqname = seqname,log_file=log_file)
        if gwexons != []:
            return gwexons
    return [[start,end,pseudo]]

#exon_finder(4179626,4179670,"-",1,40,40,0,0,cra.genome_sequence,'Chromosome3') #debug>>> exon_finder(4181323,4181442,"-",1,40,40,2,2,cra.genome_sequence,'Chromosome3') #debug

def exontuples2gff(exontuple_list,strand,feature_name,locus_name):
    outlines = []
    problems = []
    is_pseudo = False
    for exontuple in exontuple_list:
        outlines.append('\t'.join([locus_name,'ABCENTH','CDS',str(min(exontuple[:2])),str(max(exontuple[:2])),'.',strand,
                                  '.','gene_id ' + feature_name + ';transcript_id ' + feature_name + '-RA']))
        if exontuple[2]:
            is_pseudo = True
            problems.extend(list(exontuple[2]))
    if is_pseudo:
        problems = "".join(list(set(problems)))
        return '\n'.join(outlines).replace(';',problems + ';').replace('-RA',problems + '-RA')
    else:
        return '\n'.join(outlines)



#Defines a few operations that will need to be executed in different parts of the following code. NOT independent functions!
def recover_missing_exon(strand,cluster,exon_numbers,includes_start,includes_stop, search_coords, log_file,args,exon_info_dict,
                        target_genome,locus,tstart,last_strand,working_annotation,full_pseudoexon_search ):
    log_file.write('Exons skipped, attempting to recover.\n')
    # V- I think I fixed it
    #needs to be fixed to adjust phases in case not all exons are found
    found_coords_list = []
    adjust_next_phase = 0
    for exon_number in exon_numbers:
        if includes_start and exon_number == exon_numbers[0]:
            is_start,is_stop = True,False
        elif includes_stop and exon_number == exon_numbers[-1]:
            is_start,is_stop = False,True
        else:
            is_start,is_stop = False,False
        if "query_exon_info_table" in args:
            missing_exon_info = exon_info_dict[str(cluster) + ':' + str(exon_number)]
            missing_coords = orf_finder(target_genome.genome_sequence[locus],int(missing_exon_info[3]),
                                                 int(missing_exon_info[4]),strand, int(missing_exon_info[5]),5, search_coords,
                                                 is_stop = is_stop,is_start = is_start,hmm_profile = missing_exon_info[6],
                                                 evalue = args.orf_finder_E, genewise_on_fail = full_pseudoexon_search)
        log_file.write(locus + '\t' + str(tstart) + '\t' + str(cluster) + '\t' + str(exon_number) + '\t' + str(search_coords) + '\t' + str(missing_coords) + '\n')
        if missing_coords != []:
            if found_coords_list != []:
                if ((found_coords_list[-1][1] > missing_coords[0][0] and last_strand == "+") or
                (found_coords_list[-1][0] < missing_coords[-1][1] and last_strand == "-")):
                    adjust_next_phase = ( int(missing_exon_info[3]) - int(missing_exon_info[4]) ) % 3
                    continue
            if last_strand == "+":
                missing_coords[0][0] += adjust_next_phase
            elif last_strand == "-":
                missing_coords[-1][1] -=  adjust_next_phase
            found_coords_list.extend(missing_coords)
            adjust_next_phase = 0
        else:
            adjust_next_phase = (adjust_next_phase + int(missing_exon_info[3]) - int(missing_exon_info[4]) ) % 3
    if adjust_next_phase != 0 and found_coords_list != []:
        if last_strand == "+":
            found_coords_list[-1][1] -= adjust_next_phase % 3
        elif last_strand == "-":
            found_coords_list[-1][0] += adjust_next_phase % 3
    if strand == "-":
        found_coords_list.reverse()
    if found_coords_list != []:
        working_annotation.extend(found_coords_list)
        if len(found_coords_list) == len(exon_numbers):
            return (working_annotation, True)
        else:
            return (working_annotation,"I")
    else:
        return (working_annotation, False)


def last_annotation_almost_complete(last_strand,working_name,last_exon_num,last_tend,last_tstart,tstart,
                                    tend,log_file,args,exon_info_dict,target_genome,locus,working_annotation,
                                    last_startphase,last_num_exons,full_pseudoexon_search):
    if last_strand == '-':
        is_start, is_stop, if_missing = True,False,"N"
        working_annotation, recovered = recover_missing_exon(last_strand,working_name[0].split('coord')[0],range(last_exon_num),
                                        is_start,is_stop,[max((last_tend,last_tstart)),min(max((last_tend,last_tstart)) + 2000,
                                        min((tstart,tend)))],log_file,args,exon_info_dict,target_genome,locus,tstart,last_strand,
                                        working_annotation,full_pseudoexon_search)
    elif last_strand == '+':
        is_start,is_stop, if_missing = False,True,"C"
        working_annotation, recovered = recover_missing_exon(last_strand,working_name[0].split('coord')[0],range(last_exon_num + 1, last_num_exons),
                                        is_start,is_stop,[max((last_tend,last_tstart)),min(max((last_tend,last_tstart)) + 2000,
                                        min((tstart,tend)))],log_file,args,exon_info_dict,target_genome,locus,tstart,last_strand,
                                        working_annotation,full_pseudoexon_search)
    if not recovered:
        if last_strand == '-' and last_startphase != 0:
            working_annotation[-1][1] = working_annotation[-1][1] - ((3 - last_startphase) % 3)
        if working_annotation[-1][2]:
            working_annotation[-1][2] = working_annotation[-1][2] + if_missing
        else:
            working_annotation[-1][2] = if_missing
    log_file.write('writing almost complete annotation "' + working_name[0] + '". \n' )
    print(exontuples2gff(working_annotation, last_strand,working_name[0],locus))

def main():
    parser = argparse.ArgumentParser(description='ABCENTH, Annotation Based on Conserved Exons Noticed Through Homology. Takes a table of \
                                    exon alignments with phase and length information and tries to build good annotations. \
                                    Currently assumes that hits DO NOT contain stop codons.')

    parser.add_argument('--table', help = 'Table of exon alignments. Expected to be in tblastn standard tab delimited format \
                        (-outfmt 6; query,target,blah,blah,blah,blah,qstart,qstop,tstart,tstop,evalue,score) with additional \
                        fields for query exon start phase, query exon end phase, query exon peptide length, query exon number (from 0), and highest query exon number \
                        (e.g. if their are 3 exons, they would be number 0, 1, and 2, so highest query exon number would be 2). \
                        Assumed to contain only non-overlapping hits (use HitTabFilter.py to find best non-overlapping hits)')

    parser.add_argument('--genome', help = 'Genome for which to generate annotations (needed to ID splice sites)')

    parser.add_argument('--sep', default = 'tab', help = 'Hit table field delimiter, accepts "comma" and "tab". Default= "tab"')

    parser.add_argument('--maxintron', type = int, default = 100000, help = 'Maximum lenth of introns. Default = 100000')

    parser.add_argument('--query_exon_info_table', help = 'Information about each query exon (for finding missing exons)')

    parser.add_argument('--full_pseudoexon_search', default = True, type = isyes, 
                        help = 'Experimental (in progress): For truncated exons, performs genewise search \
                        against exon hmm to fill out exon in case of frameshift or stop codon')

    parser.add_argument('--orf_finder_E', type = float, default = 0.05, 
                        help = 'minimum hmmer E value to keep exons found through OR/phase matching')

    parser.add_argument('--log', default = None, help = "file to write a log for debugging (optional)")

    args = parser.parse_args()

    if args.log:
        log_file = open(args.log,'w')
    else:
        log_file = open(os.devnull, 'w')

    if args.sep == "comma":
        delimiter = ','
    elif args.sep == "tab":
        delimiter = "\t"

    if "query_exon_info_table" in args:
        exon_info_dict = {}
        for line in open(args.query_exon_info_table):
            fields = line.replace('\n','').replace('\r','').split('\t')
            exon_id = fields[0] + ':' + fields[1]
            exon_info_dict[exon_id] = fields



    #Builds nested dictionaries for sorting by locus and co-ordinates
    locus_dict = {}

    for line in open(args.table):
        fields = line.replace('\r','').replace('\n','').split(delimiter)
        if not fields[1] in locus_dict:
            locus_dict[fields[1]] = {}
        locus_dict[fields[1]][int(fields[8])] = fields

    #Imports genome to evaluate splice sites and reading frames
    target_genome = genome.Genome(args.genome,truncate_names = True)

    #Here we go
    annotation_list = []
    for locus in locus_dict.keys():
        log_file.write('starting annotation of locus: ' + locus + '\n')
        working_annotation = []
        working_name = ["oops",0]
        #sets up variables needed for determining whether or not to start a new annotation
        last_from_start,last_from_end,last_strand,last_exon_num,last_num_exons,last_startphase,last_endphase,last_tstart,last_tend = 0,None,None,None,None,None,None,0,0
        #sorts coords
        sorted_coords = sorted(locus_dict[locus].keys())
        #goes through all hits on locus
        for start_coord_index in range(len(locus_dict[locus].keys()) + 1):
            #assigns values if not the first iteration (which is only used to get values for next iteration)
            if not start_coord_index == 0:
                qstart,qend,tstart,tend,startphase,endphase,qlen,exon_num,num_exons,strand,from_start,from_end = \
                    nextqstart,nextqstop,nexttstart,nexttstop,nextstartphase,nextendphase,nextqlen,nextexon_num,\
                    nextnum_exons,nextstrand,nextfrom_start,nextfrom_end
                hit = locus_dict[locus][sorted_coords[start_coord_index - 1]]
            #gets values for next iteration and to compare the comming hit to the current one and the previous one
            if start_coord_index < len(locus_dict[locus].keys()):
                next_start_coord = sorted_coords[start_coord_index]
                nexthit = locus_dict[locus][next_start_coord]
                nextqstart,nextqstop,nexttstart,nexttstop = int(nexthit[6]),int(nexthit[7]),int(nexthit[8]),int(nexthit[9])
                if nexttstart > nexttstop:
                    nextstrand = '-'
                else:
                    nextstrand = '+'
                nextstartphase,nextendphase,nextqlen,nextexon_num,nextnum_exons = int(nexthit[12]),\
                int(nexthit[13]), int(nexthit[14]) ,int(nexthit[15]),int(nexthit[16])
                if nextstrand == '-':
                    nextfrom_start = nextnum_exons - nextexon_num
                    nextfrom_end = nextexon_num
                else:
                    nextfrom_start = nextexon_num
                    nextfrom_end = nextnum_exons - nextexon_num
            #Has program move on if this is the first iteration (which is only used to get values for next iteration)
            if start_coord_index == 0:
                continue
            #Determines if exon is good and whether or not to start a new annotation

            if from_start == 0 and (nextfrom_start == 1 or working_annotation == [] or last_from_start > nextfrom_start or last_strand != nextstrand == strand):
                #found first/last exon! Starting new annotation
                if len(working_annotation) != 0:
                    #last annotation had something! Bulding an incomplete gene model for it (this can be scrapped later, after all)
                    last_annotation_almost_complete(last_strand,working_name,last_exon_num,last_tend,last_tstart,tstart,
                                    tend,log_file,args,exon_info_dict,target_genome,locus,working_annotation,
                                    last_startphase,last_num_exons,args.full_pseudoexon_search)
                working_name = [hit[0] +"coords" + locus + "-" + str(tstart),float(hit[11])]
                log_file.write('found start/end of gene, beggining new annotation near ' + str(tstart) + ' with working name "' + working_name[0] + '\n')
                if strand == "+":
                    working_annotation = exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,qend,
                        qlen,startphase,endphase,target_genome.genome_sequence,locus,is_start = True,
                        cluster = working_name[0].split('coord')[0],exon_number = exon_num, log_file = log_file,
                        full_pseudoexon_search = args.full_pseudoexon_search, exon_info_dict = exon_info_dict)
                elif strand == '-':
                    working_annotation = exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,qend,
                        qlen,startphase,endphase,target_genome.genome_sequence,locus,is_stop = True,
                        cluster = working_name[0].split('coord')[0],exon_number = exon_num, log_file = log_file,
                        full_pseudoexon_search = args.full_pseudoexon_search, exon_info_dict = exon_info_dict)
        
            elif (from_start == 1 and (working_annotation == [] or tstart - last_tstart > args.maxintron or
                                    (nextstrand != last_strand == strand or last_from_start > nextfrom_start > from_start or nextnum_exons == num_exons != last_num_exons))
                ) or last_from_start > nextfrom_start > from_start or last_strand != nextstrand == strand or tstart - last_tstart > args.maxintron:
                #looks like we've hit a new gene even though we didn't get the first exon. Starting new annotation
                if len(working_annotation) != 0:
                    #last annotation had something! Bulding an incomplete gene model for it (this can be scrapped later, after all)
                    last_annotation_almost_complete(last_strand,working_name,last_exon_num,last_tend,last_tstart,tstart,
                                    tend,log_file,args,exon_info_dict,target_genome,locus,working_annotation,
                                    last_startphase,last_num_exons,args.full_pseudoexon_search)
                working_name = [hit[0] +"coords" + locus + "-" + str(tstart),float(hit[11])]
                working_annotation = []
                log_file.write('found almost start/end of gene, beggining new annotation near ' + str(tstart) + ' with working name "' + working_name[0] + '\n')
                #OK, now trying to recover first/last exon
                if strand:
                    if strand == '-':
                        working_annotation, recovered = recover_missing_exon(strand,working_name[0].split('coord')[0],range(exon_num + 1,num_exons + 1),
                                                        False,True, [max((min((tend,tstart)) - 2000,max((last_tend,last_tstart)))),
                                                        min((tend,tstart))],log_file, args,exon_info_dict,target_genome,locus,tstart,last_strand,
                                                        working_annotation,args.full_pseudoexon_search)
                    elif strand == '+':
                        working_annotation, recovered = recover_missing_exon(strand,working_name[0].split('coord')[0],range(exon_num),
                                                        True,False, [max((min((tend,tstart)) - 2000,max((last_tend,last_tstart)))),
                                                        min((tend,tstart))],log_file,args,exon_info_dict,target_genome,locus,tstart,last_strand,
                                                        working_annotation,args.full_pseudoexon_search)
                else:
                    recovered = False
                #Now adding found exon
                if strand == "+" and working_annotation == []:
                    working_annotation.extend(exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,
                        qend,qlen,0,endphase,target_genome.genome_sequence,locus,is_start = True,nevermind_atg = True,
                        cluster = working_name[0].split('coord')[0],exon_number = exon_num, log_file = log_file,
                        full_pseudoexon_search = args.full_pseudoexon_search, exon_info_dict = exon_info_dict))
                    if_fail = "N"
                else:
                    working_annotation.extend(exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,
                        qend,qlen,startphase,endphase,target_genome.genome_sequence,locus,
                        cluster = working_name[0].split('coord')[0],exon_number = exon_num, log_file = log_file,
                        full_pseudoexon_search = args.full_pseudoexon_search, exon_info_dict = exon_info_dict))
                    if_fail = "C"
                if len(working_annotation) < 2:
                    if working_annotation[0][2]:
                        working_annotation[0][2] = working_annotation[0][2] + if_fail #Flags as incomplete because first/last exon was missed
                    else:
                        working_annotation[0][2] = if_fail
            ##Shouldn't need this since I'm starting a new annotation anytime the intron was too long
            # elif tstart - last_tstart > args.maxintron:
            #     #Intron was too long! Clearing annotation (or writing if almost done)
            #     if last_from_end == 1 and working_annotation != []:
            #         last_annotation_almost_complete(last_strand,working_name,last_exon_num,last_tend,last_tstart,tstart,
            #                        tend,log_file,args,exon_info_dict,target_genome,locus,working_annotation,
            #                        last_startphase,last_num_exons,args.full_pseudoexon_search)
            #     working_annotation = []
            #     working_name = ["oops",0]
            elif num_exons == last_num_exons and working_annotation != [] and strand == last_strand and from_start > last_from_start:
                #Found the next exon in the progression!
                if float(hit[11]) > working_name[1]:
                    working_name == [hit[0] +"coords" + locus + "-" + str(tstart),float(hit[11])]
                if from_start > last_from_start + 1:
                    #One or more exons were skipped! Attempting to recover
                    working_annotation, recovered = recover_missing_exon(strand,working_name[0].split('coord')[0],range(min([last_exon_num + 1,exon_num + 1]),
                                                                        max([last_exon_num ,exon_num])),False,False,[max((last_tend,last_tstart)),
                                                                        min((tstart,tend))],log_file,args,exon_info_dict,target_genome,locus,tstart,last_strand,
                                                                        working_annotation,args.full_pseudoexon_search)
                if from_end != 0:
                    working_annotation.extend(exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,qend,
                        qlen,startphase,endphase,target_genome.genome_sequence,locus,
                        cluster = working_name[0].split('coord')[0],exon_number = exon_num, log_file = log_file,
                        full_pseudoexon_search = args.full_pseudoexon_search, exon_info_dict = exon_info_dict))
                    if from_start > last_from_start + 1 and recovered != True:
                        #need to adjust exon start to account for mis-matched phase
                        if working_annotation[-1][2]:
                            working_annotation[-1][2] = working_annotation[-1][2] + "I"
                        else:
                            working_annotation[-1][2] = "I"
                        if not recovered:
                            if strand == "+":
                                exon_start_adjust = (last_endphase - startphase) % 3
                                working_annotation[-1][0] = working_annotation[-1][0] + exon_start_adjust
                                log_file.write('adjusting phase for missing exon, last end-phase = ' + str(last_endphase) + ', next start-phase =' +
                                            str(startphase) + ', exon-adjust = ' + str(exon_start_adjust) + '\n')
                            elif strand == "-":
                                exon_start_adjust = (endphase - last_startphase) % 3
                                working_annotation[-2][1] = working_annotation[-2][1] - exon_start_adjust
                                log_file.write('adjusting phase for missing exon, last end-phase = ' + str(last_endphase) + ', next start-phase =' +
                                            str(startphase) + ', exon-adjust = ' + str(exon_start_adjust) + '\n')
                else:
                    if strand == "+":
                        working_annotation.extend(exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,qend,
                            qlen,startphase,endphase,target_genome.genome_sequence,locus, is_stop = True,
                            cluster = working_name[0].split('coord')[0],exon_number = exon_num, log_file = log_file,
                            full_pseudoexon_search = args.full_pseudoexon_search, exon_info_dict = exon_info_dict))
                    elif strand == '-':
                        working_annotation.extend(exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,qend,
                            qlen,startphase,endphase,target_genome.genome_sequence,locus, is_start = True,
                            cluster = working_name[0].split('coord')[0],exon_number = exon_num, log_file = log_file,
                            full_pseudoexon_search = args.full_pseudoexon_search, exon_info_dict = exon_info_dict))
                    if from_start > last_from_start + 1 and not recovered == True:
                        #need to adjust exon start to account for mis-matched phase
                        if working_annotation[-1][2]:
                            working_annotation[-1][2] = working_annotation[-1][2] + "I"
                        else:
                            working_annotation[-1][2] = "I"
                        if not recovered:
                            if strand == "+":
                                exon_start_adjust = (last_endphase - startphase) % 3
                                working_annotation[-1][0] = working_annotation[-1][0] + exon_start_adjust
                                log_file.write('adjusting phase for missing exon, last end-phase = ' + str(last_endphase) + ', next start-phase =' +
                                            str(startphase) + ', exon-adjust = ' + str(exon_start_adjust) + '\n')
                            elif strand == "-":
                                exon_start_adjust = (endphase - last_startphase) % 3
                                working_annotation[-2][1] = working_annotation[-2][1] - exon_start_adjust
                                log_file.write('adjusting phase for missing exon, last end-phase = ' + str(last_endphase) + ', next start-phase =' +
                                            str(startphase) + ', exon-adjust = ' + str(exon_start_adjust) + '\n')
                    print(exontuples2gff(working_annotation, strand,working_name[0],locus))
                    log_file.write('writing complete annotation "' + working_name[0] + '"\n')
                    working_annotation = []
            else:
                continue
            last_from_start,last_from_end,last_strand,last_exon_num,last_num_exons,last_startphase,last_endphase,last_tstart,last_tend = \
                from_start,from_end,strand,exon_num,num_exons,startphase,endphase,tstart,tend
        if len(working_annotation) != 0:
            #last annotation was almost there! Bulding an incomplete gene model for it (this can be scrapped later, after all)
            last_annotation_almost_complete(last_strand,working_name,last_exon_num,last_tend,last_tstart,tstart,
                                    tend,log_file,args,exon_info_dict,target_genome,locus,working_annotation,
                                    last_startphase,last_num_exons,args.full_pseudoexon_search)

if __name__ == "__main__":
    main()