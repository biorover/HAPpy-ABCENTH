#!/usr/bin/env python

#Wraps hmmsearch to give the user a crude ability to conduct an HMMER search of a protein query
#against a nucleotide database.

#One day, one glorious day, the HMMER developers will write a tool to do this correctly. If you're frustrated
#with this script, join with me to bug them to add this feature themselves. I hear rumors that HMMER4 will be out
#soon- would it be too much to hope for that it will contain this function?

import subprocess
from HAPpy import genome_fork as genome
import argparse
import os
import tempfile
import shlex
import sys


def main():
    parser = argparse.ArgumentParser(description='thammerin: Pronounced Tee hammerin\', so pronounced identically \
                                    to tHMMERn but does NOT infinge on the copyrights of HHMI (HMMER) or \
                                    NCBI (tblastn). \nWraps hmmsearch to give the user a crude ability to conduct \
                                    a HMMER search of a protein query against a nucleotide database. \n\nProgram\
                                    dependencies:python (obviously), hmmer suite v3, and the genome_fork library from HAP.py. \n\
                                    Nota bene: for a 200 Mbp genome, building a database of translated sequences takes \
                                    ca. 5 minutes and running hmmsearch against this database takes ca. 5 seconds. Therefore \
                                    if you are going to run a search on the same genome multiple times it will GREATLY improve \
                                    perfomance if you save the database from the first run with the "--frames_out" option \
                                    and use it in subsequent runs with the "--frames_in" option.')

    parser.add_argument('-n','--nucleotide_seqs',dest='target_nucdb',
                        help = 'Nucleotide sequences to be searched (in fasta format)')
    parser.add_argument('-p','--protein_hmm', dest = 'hmm_file',
                        help = 'Query protein hmms (created using the hmmbuild program)')
    parser.add_argument('-e','--iEvalue_threshold', dest = 'iEval_thresh', default = float('inf'), type=float,
                        help = 'threshold "i-Evalue" for inclusion in results')
    parser.add_argument('-s','--min_orf_size', dest = 'min_orf_size', default = 10,
                        help = 'minimum orf to search with HMM profile')
    parser.add_argument('-f','--hmmsearch_filepath',dest = 'hmmsearch_filepath', default = 'hmmsearch',
                        help = 'optional path to hmmsearch program if not in your PATH variable')
    # parser.add_argument('--out_format', dest = 'out_format', default = 'standard',
    #                     help = 'output format, currently defunct, may eventually accept any combination of hmm name, target name, percent identity,\
    #                     alignment length, mismatch number, null field (for compatibility with blast standard output), hmm start,\
    #                     hmm end, target start, target end, hsp evalue, hmmsearch score, standard (all of the preceeding fields in\
    #                     that order), number of idenitical positions, number of similar positions, number of gaps, \
    #                     percentage of similar positions, any specific string (in quotes), target frame, hmm match sequence,\
    #                     target match sequeence, hmm length, and target length (whole seq, not just frame)')
    parser.add_argument('--frames_out' , default = False, help = 'File to save the translated frames \
                        if you want to keep them (e.g. to use for a later run)')
    parser.add_argument('--frames_in', default = False, help = 'File with frames from a previous run (will not re-make frames)')
    parser.add_argument('--hmmsearch_log', default = False, help = 'Log file to write the full hmmsearch results to (optional)')
    parser.add_argument('--add_lengths', default = True,type = bool, help = 'Whether to add qlen and slen to stdout (default = True)')

    args=parser.parse_args()

    #sets up default program paths, overwritten by any program paths passes with -f or --program_filepaths
    hmmsearch = args.hmmsearch_filepath
    target_nucdb = genome.Genome(args.target_nucdb,truncate_names = True)

    #Gets ORFs from the genome and hmmers them

    if args.frames_in:
        frames_file = args.frames_in
    else:
        if args.frames_out:
            frame_fasta = open(args.frames_out,'w')
        else:
            frame_fasta = tempfile.NamedTemporaryFile('w')
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
                    if len(orf) > args.min_orf_size:
                        fasta_list.append('>' + seq_id + "_frameOffset-" + str(frame_offset) + "_orfStart-" + str(orf_start) +
                                        '\n' + orf)
            frame_fasta.write('\n'.join(fasta_list) + '\n')
        frame_fasta.flush()
        frames_file = frame_fasta.name

    output_tempfile = tempfile.NamedTemporaryFile()
    output_logfile = tempfile.NamedTemporaryFile()
    subprocess.call(shlex.split(hmmsearch + " --domtblout " + output_tempfile.name + " " + args.hmm_file + " " + frames_file ), stdout = output_logfile)
    hmmresults = open(output_tempfile.name)
    for line in hmmresults:
        if line[0] != "#":
            fields = line.split()
            if float(fields[12]) < args.iEval_thresh:
                tname = fields[0]
                frame_offset = int(tname.split('_frameOffset-')[1].split('_orfStart-')[0])
                orf_start = int(tname.split('_orfStart-')[1])
                if frame_offset < 3:
                    strand = "+"
                    start = orf_start * 3 + frame_offset + int(fields[17]) * 3 - 2
                    stop = orf_start * 3 + frame_offset + int(fields[18]) * 3
                else:
                    strand = "-"
                    start = frame_offset - 3 * orf_start - 3 * int(fields[17]) + 3
                    stop = frame_offset - 3 * orf_start - 3 * int(fields[18]) + 1
                if args.add_lengths:
                    fields[13] += "\t" + fields[5] + '\t' + fields[2]
                sys.stdout.write("\t".join([fields[3],tname.split('_frameOffset-')[0],fields[21],'.','.',strand,
                                fields[15],fields[16],str(start),str(stop),fields[12],fields[13]]) + '\n')
    if args.hmmsearch_log:
        output_logfile.flush()
        hmmsearch_log = open(args.hmmsearch_log,'w')
        for line in open(output_logfile):
            hmmsearch_log.write(line)
        hmmsearch_log.close()
    try:
        output_tempfile.close()
        output_logfile.close()
        frame_fasta.close()
    except:
        pass

if __name__ == "__main__":
    main()