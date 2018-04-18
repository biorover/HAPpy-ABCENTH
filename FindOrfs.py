#!/usr/bin/env python

import argparse
import genome_fork as genome
import re
import subprocess
import tempfile

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Finds orfs between splice sites (or start and stop sites) which match a given \
                                     starting and ending intron phase')
    
    parser.add_argument('--sequence', help = 'sequence, raw or in fasta format')
    parser.add_argument('--startphase', type = int, help = 'phase of preceding intron')
    parser.add_argument('--stopphase', type = int, help = 'phase of following intron')
    parser.add_argument('--expected_aa_length', type = int, help = 'expected length (in amino acids) of exon product')
    parser.add_argument('--length_variance', help = 'permitted variance in exon length (in amino acids)')
    parser.add_argument('--strand', help = 'strand of expected feature')
    parser.add_argument('--search_coords', default = None, help = 'subset of sequence within which to search')
    parser.add_argument('--is_start', default = False, type = bool, help = 'Is this exon the first exon- i.e. should it start with atg instead of a splice site')
    parser.add_argument('--is_stop', default = False, type = bool, help = 'Is this exon the last exon- i.e. should it start with a stop codon instead of a splice site')
    parser.add_argument('--hmm_profile', default = None, help = "If provided, FindOrfs will search orfs against the provided\
                        hmm profile and return the best hit")
    
    args = parser.parse_args()

def orf_finder(sequence,startphase,stopphase,strand,expected_aa_len,length_variance,
               search_coords = None, is_start = False, is_stop = False, hmm_profile = None):
    """Finds all putative exons matching a given expected length and intron phase profile and \
    containing an open reading frame"""
    search_seq = sequence
    if search_coords:
        search_seq = sequence[search_coords[0]:search_coords[1]]
    if strand == '+':
        search_seq = genome.Sequence(search_seq)
    elif strand == '-':
        search_seq = genome.Sequence(search_seq).reverse_compliment()
    match_len = expected_aa_len * 3 + ((3 - startphase) % 3) + stopphase
    phase_matches = []
    startstop = ['AG','G[TC]']
    if is_start:
        startstop[0] = 'ATG'
        match_len = match_len - 3
    if is_stop:
        startstop[1] = "T(AG|GA|AA)"
        match_len = match_len - 3
    for variance in range(length_variance + 1):
        for direction in (1,-1):
            if variance == 0 and direction == -1: # prevents double returns of 0 variance
                continue
            for match in re.finditer(startstop[0] + '.{' + str(match_len + variance * 3 * direction) +
                                     '}' + startstop[1], search_seq):
                phase_matches.append([match.start(), match.start() + match_len + variance * 3 *
                                      direction + len(startstop[0]) + len(startstop[1])])
    exon_coords = []
    for match_coords in phase_matches:
        start = match_coords[0] + len(startstop[0]) 
        stop = match_coords[1] - len(startstop[1])
        if stop - (start + ((3 - startphase) % 3)) > 2:
            if not "*" in genome.Sequence(search_seq[start + ((3 - startphase) % 3) : stop ]).translate():
                #NB: returns 1-based coords consistent with gff and blast coords
                exon_coords.append([start + 1,stop])
                if is_start:
                    exon_coords[-1][0] = exon_coords[-1][0] - 3
                if is_stop:
                    exon_coords[-1][1] = exon_coords[-1][1] + 3
                if strand == '-':
                    exon_coords[-1] = [len(search_seq) - exon_coords[-1][1] + 1,len(search_seq) - exon_coords[-1][0] + 1]
                if search_coords:
                    exon_coords[-1] = [exon_coords[-1][0] + search_coords[0],exon_coords[-1][1] + search_coords[0]]     
    if hmm_profile and len(exon_coords) > 0:
        exon_coords = hmmsearch(hmm_profile,exon_coords,sequence,strand,startphase)
    return exon_coords

def hmmsearch(hmm_profile,exon_coords,sequence,strand,startphase):
    orf_file = tempfile.NamedTemporaryFile('w')
    for coords_index in range(len(exon_coords)):
        nuc_seq = genome.Sequence(sequence[exon_coords[coords_index][0] - 1:exon_coords[coords_index][1]])
        if strand == '-':
            nuc_seq = nuc_seq.reverse_compliment()
        pep_seq = genome.Sequence(nuc_seq[((3 - startphase) % 3):]).translate()
        orf_file.write(">coords" + str(coords_index) + '\n' +
                       pep_seq + '\n' )
    orf_file.flush()
    hmmout = subprocess.check_output('hmmsearch --max -E 0.5 ' + hmm_profile + ' ' + orf_file.name,
                                     shell = True).split('\n')
    found_hit = False
    for line in hmmout:
        if line[:2] == ">>":
            found_hit = True
            return exon_coords[int(line[9:].replace('\r',''))]
            break
    if not found_hit:
        return []


def main():
    sequence = ""
    seqcounter = 0
    for line in open(args.sequence):
        if not line[0] == ">":
            sequence = sequence + line.replace('\r','').replace('\n','')
        else:
            seqcounter += 1
            if seqcounter > 1:
                print "multiple seqs in fasta- not yet supported"
                break    
    exon_coords = orf_finder(sequence,args.startphase, args.stopphase, args.strand,
                             args.expected_aa_length, args.length_variance, args.search_coords,
                             args.is_start, args.is_stop)
    for coords in exon_coords:
        print coords

if __name__ == "__main__":
    main()