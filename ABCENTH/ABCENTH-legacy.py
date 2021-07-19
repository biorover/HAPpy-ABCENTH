#!/usr/bin/env python

#   Features to be added later:
#proper splice site modeler (perhaps optimizing splice site score vs. difference from expected length), maybe MAYBE taking
#into account non-conserved exon phase.

#   Known issues:
#Will return annotations with arbitrary missing exons, but can't handle two exons missing in a row

import argparse
import genome_fork as genome
from FindOrfs import orf_finder


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

parser.add_argument('--log', default = None, help = "file to write a log for debugging (optional)")

args = parser.parse_args()


if args.log:
    log_file = open(args.log,'w')

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


#The heart of gold of this program- the function to define an exon based on the coords and target phase
def exon_finder(tstart,tend,strand,qstart,qend,qlen,qstartphase,qendphase,seqdict,seqname,
                max_offset = 30, is_start = False, is_stop = False, nevermind_atg = False):
    """"finds exons with ORFs based on the requested parameters. Note that it is expected that tstart < tend, \
    even for minus strand features, so these might need to be reversed if coming from say tblastn. \
    Also coords are expected as 1-based (as output from blast), and are \
    converted internally to 0 based. It is also expected that the hit itself doesn't contain stop codons."""
    start = None
    end = None
    pseudo = False
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
    return [start,end,pseudo]

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


#Builds nested dictionaries for sorting by locus and co-ordinates
locus_dict = {}

for line in open(args.table):
    fields = line.replace('\r','').replace('\n','').split(delimiter)
    if not fields[1] in locus_dict:
        locus_dict[fields[1]] = {}
    locus_dict[fields[1]][int(fields[8])] = fields

#Imports genome to evaluate splice sites and reading frames
target_genome = genome.Genome(args.genome)

#Defines a few operations that will need to be executed in different parts of the following code. NOT independent functions!
def recover_missing_exon(strand,cluster,exon_number,is_start,is_stop, search_coords):
    if "query_exon_info_table" in args:
        missing_exon_info = exon_info_dict[str(cluster) + ':' + str(exon_number)]
        missing_coords = orf_finder(target_genome.genome_sequence[locus],int(missing_exon_info[3]),
                                             int(missing_exon_info[4]),strand, int(missing_exon_info[5]),5, search_coords,
                                             is_stop = is_stop,is_start = is_start,hmm_profile = missing_exon_info[6])
    if args.log:
        log_file.write(locus + '\t' + str(tstart) + '\t' + str(cluster) + '\t' + str(exon_number) + '\t' + str(search_coords) + '\t' + str(missing_coords) + '\n')
    if missing_coords != []:
        working_annotation.append(missing_coords + [False])
        return True
    else:
        return False


def last_annotation_almost_complete():
    if last_strand == '-':
        is_start, is_stop, if_missing = True,False,"N"
    elif last_strand == '+':
        is_start,is_stop, if_missing = False,True,"C"
    recovered = recover_missing_exon(last_strand,working_name[0].split('coord')[0],eval(str(last_exon_num) + last_strand + "1"),is_start,is_stop,[max((last_tend,last_tstart)),min(max((last_tend,last_tstart)) + 2000,min((tstart,tend)))])
    if not recovered:
        if last_strand == '-' and last_startphase != 0:
            working_annotation[-1][1] = working_annotation[-1][1] - ((3 - last_startphase) % 3)
        if working_annotation[-1][2]:
            working_annotation[-1][2] = working_annotation[-1][2] + if_missing
        else:
            working_annotation[-1][2] = if_missing
    print(exontuples2gff(working_annotation, last_strand,working_name[0],locus))


#Here we go
annotation_list = []
for locus in locus_dict.keys():
    working_annotation = []
    working_name = ["oops",0]
    #sets up variables needed for determining whether or not to start a new annotation
    last_from_start,last_from_end,last_strand,last_exon_num,last_num_exons,last_startphase,last_endphase,last_tstart,last_tend = None,None,None,None,None,None,None,0,0
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
        if from_start == 0 and (nextfrom_start == 1 or working_annotation == []):
            #found first/last exon! Starting new annotation
            if last_from_end == 1 and len(working_annotation) != 0:
                #last annotation was almost there! Bulding an incomplete gene model for it (this can be scrapped later, after all)
                last_annotation_almost_complete()
            working_name = [hit[0] +"coords" + locus + "-" + str(tstart),float(hit[11])]
            if strand == "+":
                working_annotation = [exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,qend,qlen,startphase,endphase,target_genome.genome_sequence,locus,is_start = True)]
            elif strand == '-':
                working_annotation = [exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,qend,qlen,startphase,endphase,target_genome.genome_sequence,locus,is_stop = True)]
        elif from_start == 1 and (working_annotation == [] or tstart - last_tstart > args.maxintron or
                                  (nextstrand != last_strand or not nextfrom_start - last_from_start in [1,2] or not nextnum_exons == last_num_exons)):
            #found second/second-to-last exon and it seems like the first/last one was missing. Starting new annotation
            if last_from_end == 1 and len(working_annotation) != 0:
                #last annotation was almost there! Bulding an incomplete gene model for it (this can be scrapped later, after all)
                last_annotation_almost_complete()
            working_name = [hit[0] +"coords" + locus + "-" + str(tstart),float(hit[11])]
            working_annotation = []
            #OK, now trying to recover first/last exon
            if strand:
                if strand == '-':
                    me_is_start, me_is_stop = False,True
                elif strand == '+':
                    me_is_start,me_is_stop =  True,False
                recovered = recover_missing_exon(strand,working_name[0].split('coord')[0],eval(str(exon_num) + strand + "-1"),me_is_start,me_is_stop, [max((min((tend,tstart)) - 2000,max((last_tend,last_tstart)))),min((tend,tstart))])
            else:
                recovered = False
            #Now adding found exon
            if strand == "+" and working_annotation == []:
                working_annotation.append(exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,qend,qlen,startphase,endphase,target_genome.genome_sequence,locus,is_start = True,nevermind_atg = True))
                if_fail = "N"
            else:
                working_annotation.append(exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,qend,qlen,startphase,endphase,target_genome.genome_sequence,locus))
                if_fail = "C"
            if len(working_annotation) < 2:
                if working_annotation[0][2]:
                    working_annotation[0][2] = working_annotation[0][2] + if_fail #Flags as pseudogene because first/last exon was missed
                else:
                    working_annotation[0][2] = if_fail
        elif tstart - last_tstart > args.maxintron:
            #Intron was too long! Clearing annotation (or writing if almost done)
            if last_from_end == 1 and working_annotation != []:
                last_annotation_almost_complete()
            working_annotation = []
            working_name = ["oops",0]
        elif num_exons == last_num_exons and working_annotation != [] and strand == last_strand and (from_start == last_from_start + 1 or (from_start == last_from_start + 2 and nextfrom_start != last_from_start + 1)):
            #Found the next exon in the progression!
            if float(hit[11]) > working_name[1]:
                working_name == [hit[0] +"coords" + locus + "-" + str(tstart),float(hit[11])]
            if from_start == last_from_start + 2:
                #An exon was skipped! Attempting to recover
                recovered = recover_missing_exon(strand,working_name[0].split('coord')[0],eval(str(exon_num) + strand + "-1"),False,False,[max((last_tend,last_tstart)),min((tstart,tend))])
            if from_end != 0:
                working_annotation.append(exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,qend,qlen,startphase,endphase,target_genome.genome_sequence,locus))
                if from_start == last_from_start + 2 and not recovered:
                    #need to adjust exon start to account for mis-matched phase
                    if working_annotation[-1][2]:
                        working_annotation[-1][2] = working_annotation[-1][2] + "I"
                    else:
                        working_annotation[-1][2] = "I"
                    if strand == "+":
                        exon_start_adjust = (last_endphase - startphase) % 3
                        working_annotation[-1][0] = working_annotation[-1][0] + exon_start_adjust
                    elif strand == "-":
                        exon_start_adjust = (endphase - last_startphase) % 3
                        working_annotation[-2][1] = working_annotation[-2][1] - exon_start_adjust
            else:
                if strand == "+":
                    working_annotation.append(exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,qend,qlen,startphase,endphase,target_genome.genome_sequence,locus, is_stop = True))
                elif strand == '-':
                    working_annotation.append(exon_finder(min((tstart,tend)),max((tstart,tend)),strand,qstart,qend,qlen,startphase,endphase,target_genome.genome_sequence,locus, is_start = True))
                if from_start == last_from_start + 2 and not recovered:
                    #need to adjust exon start to account for mis-matched phase
                    if working_annotation[-1][2]:
                        working_annotation[-1][2] = working_annotation[-1][2] + "I"
                    else:
                        working_annotation[-1][2] = "I"
                    if strand == "+":
                        exon_start_adjust = (last_endphase - startphase) % 3
                        working_annotation[-1][0] = working_annotation[-1][0] + exon_start_adjust
                    elif strand == "-":
                        exon_start_adjust = (endphase - last_startphase) % 3
                        working_annotation[-2][1] = working_annotation[-2][1] - exon_start_adjust
                print(exontuples2gff(working_annotation, strand,working_name[0],locus))
                working_annotation = []
        else:
            continue
        last_from_start,last_from_end,last_strand,last_exon_num,last_num_exons,last_startphase,last_endphase,last_tstart,last_tend = \
            from_start,from_end,strand,exon_num,num_exons,startphase,endphase,tstart,tend
    if last_from_end == 1 and len(working_annotation) != 0:
                #last annotation was almost there! Bulding an incomplete gene model for it (this can be scrapped later, after all)
                last_annotation_almost_complete()
